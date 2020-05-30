/*! \file record-full.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "constants.h"
#include "common.h"
#include "audio_logger.h"

#include <fstream>
#include <chrono>
#include <thread>

bool g_terminate = false;

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

    bool doRecord = true;
    size_t totalSize_bytes = 0;

    std::ofstream fout(argv[1], std::ios::binary);
    if (fout.good() == false) {
        fprintf(stderr, "Failed to open file '%s'\n", argv[1]);
        return -1;
    }

    AudioLogger audioLogger;
    AudioLogger::Callback cbAudio = [&](const auto & frames) {
        doRecord = true;

        for (const auto & frame : frames) {
            totalSize_bytes += sizeof(frame[0])*frame.size();
            fout.write((char *)(frame.data()), sizeof(frame[0])*frame.size());
            fout.flush();
        }

        printf("Total data saved: %g MB\n", ((float)(totalSize_bytes)/1024.0f/1024.0f));
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

    while (true) {
        if (doRecord) {
            doRecord = false;
            audioLogger.record(0.5f, 0);
        } else {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
    }

    return 0;
}
