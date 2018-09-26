/*! \file record.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "audio_logger.h"

#include <stdio.h>
#include <termios.h>
#include <unistd.h>

#include <cstdio>
#include <chrono>
#include <thread>
#include <fstream>

int main(int argc, const char ** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s output.kbd\n", argv[0]);
        return -127;
    }

    auto tStart = std::chrono::high_resolution_clock::now();
    auto tEnd = std::chrono::high_resolution_clock::now();

    size_t totalSize_bytes = 0;
    constexpr float kBufferSize_s = 0.075f;
    constexpr uint64_t kSampleRate = 24000;

    int keyPressed = -1;
    int bufferSize_frames = 2*AudioLogger::getBufferSize_frames(kSampleRate, kBufferSize_s) - 1;

    std::ofstream fout(argv[1], std::ios::binary);
    fout.write((char *)(&bufferSize_frames), sizeof(bufferSize_frames));

    AudioLogger audioLogger;
    AudioLogger::Callback cbAudio = [&](const auto & frames) {
        tEnd = std::chrono::high_resolution_clock::now();
        printf("t = %d us\n", (int) std::chrono::duration_cast<std::chrono::microseconds>(tEnd - tStart).count());

        printf("Record callback with %d frames\n", (int) frames.size());
        fout.write((char *)(&keyPressed), sizeof(keyPressed));
        for (const auto & frame : frames) {
            totalSize_bytes += sizeof(frame[0])*frame.size();
            fout.write((char *)(frame.data()), sizeof(frame[0])*frame.size());
            fout.flush();
        }
        printf("Total data saved: %g MB\n", ((float)(totalSize_bytes)/1024.0f/1024.0f));
        keyPressed = -1;
    };

    if (audioLogger.install(kSampleRate, cbAudio) == false) {
        fprintf(stderr, "Failed to install audio logger\n");
        return -1;
    }

    std::thread keyReader = std::thread([&]() {
        struct termios oldt, newt;
        tcgetattr ( STDIN_FILENO, &oldt );
        newt = oldt;
        newt.c_lflag &= ~( ICANON | ECHO );
        tcsetattr ( STDIN_FILENO, TCSANOW, &newt );
        while (true) {
            int key = getchar();
            tStart = std::chrono::high_resolution_clock::now();
            if (keyPressed == -1) {
                keyPressed = key;
                audioLogger.record(kBufferSize_s);
            }
        }
        tcsetattr ( STDIN_FILENO, TCSANOW, &oldt );
    });

    while (true) {}

    fout.close();

    return 0;
}
