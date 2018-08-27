/*! \file audio_logger.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "audio_logger.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_audio.h>

#include <mutex>
#include <algorithm>

namespace {
    void cbAudioReady(void * userData, uint8_t * stream, int32_t nbytes) {
        AudioLogger * logger = (AudioLogger *)(userData);
        logger->addFrame((AudioLogger::Sample *)(stream));
    }
}

struct AudioLogger::Data {
    SDL_AudioDeviceID deviceIdIn = 0;
    //SDL_AudioDeviceID deviceIdOut = 0;

    Callback callback;

    int32_t nFramesToRecord = 0;
    int32_t sampleSize_bytes = 4;

    int32_t bufferId = 0;
    std::array<Frame, kBufferSize_frames> buffer;

    Record record;

    std::mutex mutex;
};

AudioLogger::AudioLogger() : data_(new AudioLogger::Data()) {}

AudioLogger::~AudioLogger() {}

bool AudioLogger::install(AudioLogger::Callback callback) {
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

    captureSpec.freq = kSampleRate;
    captureSpec.format = AUDIO_F32SYS;
    captureSpec.channels = 1;
    captureSpec.samples = kSamplesPerFrame;
    captureSpec.callback = ::cbAudioReady;
    captureSpec.userdata = this;

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

    data.callback = std::move(callback);

    return true;
}

bool AudioLogger::addFrame(const Sample * stream) {
    std::lock_guard<std::mutex> lock(data_->mutex);

    auto & curFrame = data_->buffer[data_->bufferId];
    std::copy(stream, stream + kSamplesPerFrame, curFrame.data());
    if (data_->nFramesToRecord > 0) {
        data_->record.push_back(curFrame);
        if (--data_->nFramesToRecord == 0) {
            printf("Callback with %d frames\n", (int) data_->record.size());
            data_->record.clear();
        }
    }
    if (++data_->bufferId >= data_->buffer.size()) {
        data_->bufferId = 0;
    }

    return true;
}

bool AudioLogger::record() {
    std::lock_guard<std::mutex> lock(data_->mutex);

    if (data_->record.size() == 0) {
        for (size_t i = 0; i < data_->buffer.size(); ++i) {
            data_->record.push_back(data_->buffer[(data_->bufferId + i)%data_->buffer.size()]);
        }
    }

    data_->nFramesToRecord = kBufferSize_frames;

    return true;
}
