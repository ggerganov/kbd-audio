/*! \file record.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "key_logger.h"
#include "audio_logger.h"

#include <cstdio>
#include <chrono>
#include <thread>
#include <fstream>

int main(int argc, const char ** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s output.kbd\n", argv[0]);
        return -127;
    }

    size_t totalSize_bytes = 0;

    std::ofstream fout(argv[1], std::ios::binary);

    AudioLogger audioLogger;
    AudioLogger::Callback cbAudio = [&fout, &totalSize_bytes](const auto & frames) {
        printf("Record callback with %d frames\n", (int) frames.size());
        for (const auto & frame : frames) {
            totalSize_bytes += sizeof(frame[0])*frame.size();
            fout.write((char *)(frame.data()), sizeof(frame[0])*frame.size());
            fout.flush();
        }
        printf("Total data saved: %g MB\n", ((float)(totalSize_bytes)/1024.0f/1024.0f));
    };

    if (audioLogger.install(cbAudio) == false) {
        fprintf(stderr, "Failed to install audio logger\n");
        return -1;
    }

    KeyLogger keyLogger;
    KeyLogger::Callback cbKey = [&audioLogger](int key) -> void {
        //printf("%s\n", KeyLogger::codeToText(key));
        audioLogger.record();
    };

    if (keyLogger.install(cbKey) == false) {
        fprintf(stderr, "Failed to install key logger\n");
        return -2;
    }

    while (true) {}

    fout.close();

    return 0;
}
