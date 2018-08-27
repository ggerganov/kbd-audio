/*! \file main.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "key_logger.h"
#include "audio_logger.h"

#include <cstdio>
#include <chrono>
#include <thread>

int main(int , const char **) {
    AudioLogger audioLogger;
    AudioLogger::Callback cbAudio = [](const auto & frames) {
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

    while (true) {
    }

    return 0;
}
