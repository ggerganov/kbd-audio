/*! \file audio_logger.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "audio_logger.h"

struct AudioLogger::Data {
};

AudioLogger::AudioLogger() : data_(new AudioLogger::Data()) {}

bool AudioLogger::install() {
    return true;
}
