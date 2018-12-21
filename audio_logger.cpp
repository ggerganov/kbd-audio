/*! \file audio_logger.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "audio_logger.h"

#include <SDL.h>
#include <SDL_audio.h>

#include <mutex>
#include <algorithm>

namespace {
    void cbAudioReady(void * userData, uint8_t * stream, int32_t /*nbytes*/) {
        AudioLogger * logger = (AudioLogger *)(userData);
        logger->addFrame((AudioLogger::Sample *)(stream));
    }
}

struct AudioLogger::Data {
    Data() {
        for (auto & frame : buffer) {
            frame.fill(0);
        }

        record.clear();
    }

    SDL_AudioDeviceID deviceIdIn = 0;
    //SDL_AudioDeviceID deviceIdOut = 0;

    Callback callback = nullptr;

    int64_t sampleRate = kMaxSampleRate;

    int32_t nFramesToRecord = 0;
    int32_t sampleSize_bytes = 4;

    int32_t bufferId = 0;
    std::array<Frame, getBufferSize_frames(kMaxSampleRate, kMaxBufferSize_s)> buffer;

    Record record;

    std::mutex mutex;
};

AudioLogger::AudioLogger() : data_(new AudioLogger::Data()) {}

AudioLogger::~AudioLogger() {}

bool AudioLogger::install(int64_t sampleRate, AudioLogger::Callback callback, int captureId) {
    auto & data = getData();

    data.sampleRate = sampleRate;

    if (SDL_Init(SDL_INIT_AUDIO) < 0) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't initialize SDL: %s\n", SDL_GetError());
        return false;
    }

    int nDevices = SDL_GetNumAudioDevices(SDL_TRUE);
    printf("Found %d capture devices:\n", nDevices);
    for (int i = 0; i < nDevices; i++) {
        printf("    - Capture device #%d: '%s'\n", i, SDL_GetAudioDeviceName(i, SDL_TRUE));
    }

    if (captureId < 0 || captureId >= nDevices) {
        printf("Invalid capture device id selected - %d\n", captureId);
        return false;
    }

    SDL_AudioSpec captureSpec;
    SDL_zero(captureSpec);

    captureSpec.freq = data.sampleRate;
    captureSpec.format = AUDIO_F32SYS;
    captureSpec.channels = 1;
    captureSpec.samples = kSamplesPerFrame;
    captureSpec.callback = ::cbAudioReady;
    captureSpec.userdata = this;

    SDL_AudioSpec obtainedSpec;
    SDL_zero(obtainedSpec);

    printf("Attempt to open capture device %d : '%s' ...\n", captureId, SDL_GetAudioDeviceName(captureId, SDL_TRUE));
    data.deviceIdIn = SDL_OpenAudioDevice(SDL_GetAudioDeviceName(captureId, SDL_TRUE), SDL_TRUE, &captureSpec, &obtainedSpec, 0);
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

    printf("Opened capture device succesfully!\n");
    printf("    Frequency:  %d\n", obtainedSpec.freq);
    printf("    Format:     %d (%d bytes)\n", obtainedSpec.format, data.sampleSize_bytes);
    printf("    Channels:   %d\n", obtainedSpec.channels);
    printf("    Samples:    %d\n", obtainedSpec.samples);

    SDL_PauseAudioDevice(data.deviceIdIn, 0);

    data.callback = std::move(callback);

    return true;
}

bool AudioLogger::terminate() {
    auto & data = getData();

    SDL_PauseAudioDevice(data.deviceIdIn, 1);
    SDL_CloseAudioDevice(data.deviceIdIn);

    return true;
}

bool AudioLogger::addFrame(const Sample * stream) {
    auto & data = getData();

	if ((int) SDL_GetQueuedAudioSize(data.deviceIdIn) > 32*sizeof(float)*kSamplesPerFrame) {
		printf("Queue size: %d\n", SDL_GetQueuedAudioSize(data.deviceIdIn));
		SDL_ClearQueuedAudio(data.deviceIdIn);
	}

    std::lock_guard<std::mutex> lock(data.mutex);

    auto & curFrame = data.buffer[data.bufferId];
    std::copy(stream, stream + kSamplesPerFrame, curFrame.data());
    if (data.nFramesToRecord > 0) {
        data.record.push_back(curFrame);
        if (--data.nFramesToRecord == 0) {
            if (data.callback) data.callback(data.record);
            data.record.clear();
        }
    }
    if (++data.bufferId >= (int) data.buffer.size()) {
        data.bufferId = 0;
    }

    return true;
}

bool AudioLogger::record(float bufferSize_s) {
    auto & data = getData();

    if (bufferSize_s > kMaxBufferSize_s) return false;

    auto bufferSize_frames = getBufferSize_frames(data.sampleRate, bufferSize_s);

    std::lock_guard<std::mutex> lock(data.mutex);

    if (data.record.size() == 0) {
        int fStart = data.bufferId - 2 + 1;
        if (fStart < 0) fStart += data.buffer.size();
        for (size_t i = 0; i < 2 - 1; ++i) {
            data.record.push_back(data.buffer[(fStart + i)%data.buffer.size()]);
        }
    }

    data.nFramesToRecord = 2*bufferSize_frames - 2;

    return true;
}

bool AudioLogger::recordSym(float bufferSize_s) {
    auto & data = getData();

    if (bufferSize_s > kMaxBufferSize_s) return false;

    auto bufferSize_frames = getBufferSize_frames(data.sampleRate, bufferSize_s);

    std::lock_guard<std::mutex> lock(data.mutex);

    if (data.record.size() == 0) {
        int fStart = data.bufferId - bufferSize_frames + 1;
        if (fStart < 0) fStart += data.buffer.size();
        for (size_t i = 0; i < bufferSize_frames - 1; ++i) {
            data.record.push_back(data.buffer[(fStart + i)%data.buffer.size()]);
        }
    }

    data.nFramesToRecord = 2*bufferSize_frames - 2;

    return true;
}

bool AudioLogger::pause() {
    auto & data = getData();
    SDL_PauseAudioDevice(data.deviceIdIn, 1);
    return true;
}

bool AudioLogger::resume() {
    auto & data = getData();
    SDL_PauseAudioDevice(data.deviceIdIn, 0);
    return true;
}
