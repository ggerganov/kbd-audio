/*! \file main.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "key_logger.h"

#include <cstdio>

int main(int , const char **) {
    KeyLogger keyLogger;
    keyLogger.install([](int key) -> void {
        printf("%s\n", KeyLogger::codeToText(key));
    });

    return 0;
}
