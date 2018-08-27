/*! \file main.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "key_logger.h"
#include "audio_logger.h"

#include <cstdio>

int main(int , const char **) {
    AudioLogger audioLogger;

    if (audioLogger.install() == false) {
        fprintf(stderr, "Failed to install audio logger\n");
        return -1;
    }

    KeyLogger keyLogger;

    KeyLogger::Callback cbKey = [](int key) -> void {
        printf("%s\n", KeyLogger::codeToText(key));
    };

    if (keyLogger.install(cbKey) == false) {
        fprintf(stderr, "Failed to install key logger\n");
        return -2;
    }

    return 0;
}
