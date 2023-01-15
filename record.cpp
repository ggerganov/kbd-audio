/*! \file record.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "constants.h"
#include "common.h"
#include "audio-logger.h"

#include <stdio.h>
#include <termios.h>
#include <unistd.h>

#include <map>
#include <cstdio>
#include <chrono>
#include <thread>
#include <deque>
#include <fstream>

int main(int argc, char ** argv) {
    printf("Usage: %s output.kbd [-cN]\n", argv[0]);
    printf("    -cN - select capture device N\n");
    printf("    -CN - number N of capture channels N\n");
    printf("\n");

    if (argc < 2) {
        return -127;
    }

    auto argm = parseCmdArguments(argc, argv);
    int captureId = argm["c"].empty() ? 0 : std::stoi(argm["c"]);
    int nChannels = argm["C"].empty() ? 0 : std::stoi(argm["C"]);

    auto tStart = std::chrono::high_resolution_clock::now();
    auto tEnd = std::chrono::high_resolution_clock::now();

    size_t totalSize_bytes = 0;
    std::deque<int> keyPressedQueue;
    std::map<int, int> nTimes;
    printf("Recording %d frames per key press\n", kBufferSizeTrain_frames);

    std::ofstream fout(argv[1], std::ios::binary);
    fout.write((char *)(&kBufferSizeTrain_frames), sizeof(kBufferSizeTrain_frames));

    AudioLogger audioLogger;
    AudioLogger::Callback cbAudio = [&](const auto & frames) {
        tEnd = std::chrono::high_resolution_clock::now();

        int keyPressed = keyPressedQueue.front();
        keyPressedQueue.pop_front();

        fout.write((char *)(&keyPressed), sizeof(keyPressed));
        for (const auto & frame : frames) {
            totalSize_bytes += sizeof(frame[0])*frame.size();
            fout.write((char *)(frame.data()), sizeof(frame[0])*frame.size());
            fout.flush();
        }
        ++nTimes[keyPressed];

        printf("Last recorded key - %3d '%s'. Total times recorded so far - %3d. Total data saved: %g MB\n",
               keyPressed, kKeyText.at(keyPressed), nTimes[keyPressed], ((float)(totalSize_bytes)/1024.0f/1024.0f));

        keyPressed = -1;
    };

    AudioLogger::Parameters parameters;
    parameters.callback = std::move(cbAudio);
    parameters.captureId = captureId;
    parameters.nChannels = nChannels;
    parameters.sampleRate = kSampleRate;
    parameters.freqCutoff_Hz = kFreqCutoff_Hz;

    if (audioLogger.install(std::move(parameters)) == false) {
        fprintf(stderr, "Failed to install audio logger\n");
        return -1;
    }

    printf("\n");
    printf("To stop capturing press Ctrl+C\n");
    printf("On some systems (e.g. MacOS), this may not work. In this case, you have to kill the process manually\n");
    printf("\n");

    std::thread keyReader = std::thread([&]() {
        struct termios oldt, newt;
        tcgetattr ( STDIN_FILENO, &oldt );
        newt = oldt;
        newt.c_lflag &= ~( ICANON | ECHO );
        tcsetattr ( STDIN_FILENO, TCSANOW, &newt );
        while (true) {
            int key = getchar();
            tStart = std::chrono::high_resolution_clock::now();
            keyPressedQueue.push_back(key);
            if (audioLogger.record(kBufferSizeTrain_s, 2) == false) {
                fprintf(stderr, "error : failed to record\n");
            }
        }
        tcsetattr ( STDIN_FILENO, TCSANOW, &oldt );
    });

    while (true) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    keyReader.join();

    fout.close();

    return 0;
}
