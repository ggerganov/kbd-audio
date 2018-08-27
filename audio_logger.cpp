/*! \file audio_logger.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "audio_logger.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_audio.h>

namespace {
    void cbAudioReady(void * userData, uint8_t * stream, int32_t nbytes) {
        printf("Callback, %d bytes\n", nbytes);
    }
}

struct AudioLogger::Data {
    SDL_AudioDeviceID deviceIdIn = 0;
    //SDL_AudioDeviceID deviceIdOut = 0;

    int sampleSize_bytes = 4;
};

AudioLogger::AudioLogger() : data_(new AudioLogger::Data()) {}

AudioLogger::~AudioLogger() {}

bool AudioLogger::install() {
    auto & data = getData();

    if (SDL_Init(SDL_INIT_AUDIO) < 0) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't initialize SDL: %s\n", SDL_GetError());
        return false;
    }

    int nDevices = SDL_GetNumAudioDevices(SDL_TRUE);
    printf("Found %d capture devices:\n", nDevices);
    for (int i = 0; i < nDevices; i++) {
        printf("    - Capture device #%d: '%s'\n", i, SDL_GetAudioDeviceName(i, SDL_TRUE));
    }

    SDL_AudioSpec captureSpec;
    SDL_zero(captureSpec);

    captureSpec.freq = 44100;
    captureSpec.format = AUDIO_F32SYS;
    captureSpec.channels = 1;
    captureSpec.samples = 1024;
    captureSpec.callback = ::cbAudioReady;
    captureSpec.userdata = NULL;

    SDL_AudioSpec obtainedSpec;
    SDL_zero(obtainedSpec);

    data.deviceIdIn = SDL_OpenAudioDevice(NULL, SDL_TRUE, &captureSpec, &obtainedSpec, 0);
    if (!data.deviceIdIn) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't open an audio device for capture: %s!\n", SDL_GetError());
        SDL_Quit();
        return false;
    }

    data.sampleSize_bytes = 4;
    switch (obtainedSpec.format) {
        case AUDIO_U8:
        case AUDIO_S8:
            data.sampleSize_bytes = 1;
            break;
        case AUDIO_U16SYS:
        case AUDIO_S16SYS:
            data.sampleSize_bytes = 2;
            break;
        case AUDIO_S32SYS:
        case AUDIO_F32SYS:
            data.sampleSize_bytes = 4;
            break;
    }

    printf("Opened capture device %d\n", data.deviceIdIn);
    printf("    Frequency:  %d\n", obtainedSpec.freq);
    printf("    Format:     %d (%d bytes)\n", obtainedSpec.format, data.sampleSize_bytes);
    printf("    Channels:   %d\n", obtainedSpec.channels);
    printf("    Samples:    %d\n", obtainedSpec.samples);

    SDL_PauseAudioDevice(data.deviceIdIn, 0);

    return true;
}
