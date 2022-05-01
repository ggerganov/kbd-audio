/*! \file audio-logger.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "audio-logger.h"

#include <SDL.h>
#include <SDL_audio.h>

#include <mutex>
#include <atomic>
#include <algorithm>

namespace {
    void cbAudioReady(void * userData, uint8_t * stream, int32_t /*nbytes*/) {
        AudioLogger * logger = (AudioLogger *)(userData);
        logger->addFrame((AudioLogger::Sample *)(stream));
    }
}

struct AudioLogger::Data {
    Data() : isReady(false) {
        for (auto & frame : buffer) {
            frame.fill(0);
        }

        for (auto & record : records) {
            record.clear();
        }

        nFramesToRecord.fill(0);
    }

    SDL_AudioDeviceID deviceIdIn = 0;
    //SDL_AudioDeviceID deviceIdOut = 0;

    int32_t sampleSize_bytes = -1;

    int32_t bufferId = 0;
    std::array<Frame, getBufferSize_frames(kMaxSampleRate, kMaxBufferSize_s)> buffer;

    int32_t nRecords = 0;
    std::array<int32_t, kMaxRecords> nFramesToRecord;
    std::array<Record, kMaxRecords> records;

    Parameters parameters;
    TFilterCoefficients filterCoefficients;

    std::mutex mutex;
    std::atomic_bool isReady;
};

AudioLogger::AudioLogger() : data_(new AudioLogger::Data()) {}

AudioLogger::~AudioLogger() {}

bool AudioLogger::install(Parameters && parameters) {
    auto & data = getData();

    if (parameters.captureId < 0) {
        fprintf(stderr, "error : invalid captureId = %d\n", parameters.captureId);
        return false;
    }

    if (parameters.nChannels < 0) {
        fprintf(stderr, "error : invalid nChannels = %d\n", parameters.nChannels);
        return false;
    }

    if (parameters.sampleRate <= 0) {
        fprintf(stderr, "error : invalid sampleRate = %d\n", (int) parameters.sampleRate);
        return false;
    }

    if (parameters.freqCutoff_Hz <= 0) {
        fprintf(stderr, "error : invalid freqCutoff_Hz = %g\n", parameters.freqCutoff_Hz);
        return false;
    }

    static bool isInitialized = false;
    if (!isInitialized && SDL_Init(SDL_INIT_AUDIO) < 0) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't initialize SDL: %s\n", SDL_GetError());
        return false;
    }
    isInitialized = true;

    int nDevices = SDL_GetNumAudioDevices(SDL_TRUE);
    printf("Found %d capture devices:\n", nDevices);
    for (int i = 0; i < nDevices; i++) {
        printf("    - Capture device #%d: '%s'\n", i, SDL_GetAudioDeviceName(i, SDL_TRUE));
    }

    if (parameters.captureId < 0 || parameters.captureId >= nDevices) {
        fprintf(stderr, "error : invalid capture device id selected - %d\n", parameters.captureId);
        return false;
    }

    SDL_AudioSpec captureSpec;
    SDL_zero(captureSpec);

    captureSpec.freq = parameters.sampleRate;
    captureSpec.format = AUDIO_F32SYS;
    captureSpec.channels = parameters.nChannels;
    captureSpec.samples = kSamplesPerFrame;
    captureSpec.callback = ::cbAudioReady;
    captureSpec.userdata = this;

    SDL_AudioSpec obtainedSpec;
    SDL_zero(obtainedSpec);

    printf("Attempt to open capture device %d : '%s' ...\n", parameters.captureId, SDL_GetAudioDeviceName(parameters.captureId, SDL_TRUE));
    data.deviceIdIn = SDL_OpenAudioDevice(SDL_GetAudioDeviceName(parameters.captureId, SDL_TRUE), SDL_TRUE, &captureSpec, &obtainedSpec, SDL_AUDIO_ALLOW_CHANNELS_CHANGE);
    if (!data.deviceIdIn) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't open an audio device for capture: %s!\n", SDL_GetError());
        SDL_Quit();
        return false;
    }

    switch (obtainedSpec.format) {
        case AUDIO_U8:
        case AUDIO_S8:
        case AUDIO_U16SYS:
        case AUDIO_S16SYS:
        case AUDIO_S32SYS:
            {
                fprintf(stderr, "error : unsupported sample format %d\n", obtainedSpec.format);
                return false;
            }
            break;
        case AUDIO_F32SYS:
            {
                data.sampleSize_bytes = 4;
            }
            break;
    }

    printf("Opened capture device succesfully!\n");
    printf("    DeviceId:   %d\n", data.deviceIdIn);
    printf("    Frequency:  %d\n", obtainedSpec.freq);
    printf("    Format:     %d (%d bytes)\n", obtainedSpec.format, data.sampleSize_bytes);
    printf("    Channels:   %d\n", obtainedSpec.channels);
    printf("    Samples:    %d\n", obtainedSpec.samples);

    if (obtainedSpec.channels != parameters.nChannels && parameters.nChannels != 0) {
        printf("Warning: obtained number of channels (%d) does not match requested (%d)\n",
               obtainedSpec.channels, parameters.nChannels);
    }

    SDL_PauseAudioDevice(data.deviceIdIn, 0);

    parameters.nChannels = obtainedSpec.channels;

    switch (parameters.filter) {
        case EAudioFilter::None:
            {
            }
            break;
        case EAudioFilter::FirstOrderHighPass:
            {
                data.filterCoefficients = ::calculateCoefficientsFirstOrderHighPass(parameters.freqCutoff_Hz, parameters.sampleRate);
            }
            break;
        case EAudioFilter::SecondOrderButterworthHighPass:
            {
                data.filterCoefficients = ::calculateCoefficientsSecondOrderButterworthHighPass(parameters.freqCutoff_Hz, parameters.sampleRate);
            }
            break;
    };

    data.parameters = parameters;
    data.isReady = true;

    // print filter paramters
    printf("    Audio Filter: %d\n", parameters.filter);
    printf("    Cutoff frequency: %g Hz\n", parameters.freqCutoff_Hz);
    printf("Capturing audio ..\n");

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

    if (data.isReady == false) return false;

	if (SDL_GetQueuedAudioSize(data.deviceIdIn) > 32*sizeof(float)*kSamplesPerFrame) {
		printf("Queue size: %d\n", SDL_GetQueuedAudioSize(data.deviceIdIn));
		SDL_ClearQueuedAudio(data.deviceIdIn);
	}

    const float norm = 1.0/data.parameters.nChannels;

    auto & curFrame = data.buffer[data.bufferId];

    for (int i = 0; i < kSamplesPerFrame; ++i) {
        Sample x = 0;
        for (int j = 0; j < data.parameters.nChannels; ++j) {
            x += stream[i*data.parameters.nChannels + j];
        }
        curFrame[i] = x*norm;
    }

    switch (data.parameters.filter) {
        case EAudioFilter::None:
            {
            }
            break;
        case EAudioFilter::FirstOrderHighPass:
            {
                for (auto & s : curFrame) {
                    s = ::filterFirstOrderHighPass(data.filterCoefficients, s);
                }
            }
            break;
        case EAudioFilter::SecondOrderButterworthHighPass:
            {
                for (auto & s : curFrame) {
                    s = ::filterSecondOrderButterworthHighPass(data.filterCoefficients, s);
                }
            }
            break;
    }

    std::lock_guard<std::mutex> lock(data.mutex);

    for (int r = 0; r < data.nRecords; ++r) {
        auto & record = data.records[r];
        auto & nFramesToRecord = data.nFramesToRecord[r];

        if (nFramesToRecord > 0) {
            record.push_back(curFrame);
            if (--nFramesToRecord == 0) {
                if (data.parameters.callback) data.parameters.callback(record);
                record.clear();

                for (int k = r + 1; k < data.nRecords; ++k) {
                    data.records[k - 1] = std::move(data.records[k]);
                    data.nFramesToRecord[k - 1] = data.nFramesToRecord[k];
                }
                --data.nRecords;
                --r;
            }
        }
    }

    if (++data.bufferId >= (int) data.buffer.size()) {
        data.bufferId = 0;
    }

    return true;
}

bool AudioLogger::record(float bufferSize_s, int32_t nPrevFrames) {
    auto & data = getData();

    if (isValidBufferSize(bufferSize_s) == false) {
        return false;
    }

    auto bufferSize_frames = getBufferSize_frames(data.parameters.sampleRate, bufferSize_s);

    if (nPrevFrames >= bufferSize_frames) {
        fprintf(stderr, "warning : invalid previous frames in record requested - %d. max allowed is %d s\n", nPrevFrames, bufferSize_frames - 1);
        return false;
    }

    std::lock_guard<std::mutex> lock(data.mutex);

    if (data.nRecords == kMaxRecords) {
        fprintf(stderr, "warning : max number of simultaneous records %d reached\n", kMaxRecords);
        return false;
    }

    auto & record = data.records[data.nRecords];

    if (record.size() == 0) {
        int fStart = data.bufferId - nPrevFrames;
        if (fStart < 0) fStart += data.buffer.size();
        for (int i = 0; i < nPrevFrames; ++i) {
            record.push_back(data.buffer[(fStart + i)%data.buffer.size()]);
        }
    } else {
        fprintf(stderr, "warning : new record requested before last has been processed. should never happen\n");
    }

    data.nFramesToRecord[data.nRecords] = bufferSize_frames - nPrevFrames;

    ++data.nRecords;

    return true;
}

bool AudioLogger::pause() {
    auto & data = getData();
    SDL_PauseAudioDevice(data.deviceIdIn, 1);
    data.nFramesToRecord.fill(0);

    return true;
}

bool AudioLogger::resume() {
    auto & data = getData();
    SDL_PauseAudioDevice(data.deviceIdIn, 0);
    return true;
}

bool AudioLogger::isValidBufferSize(float bufferSize_s) const {
    if (bufferSize_s <= 0) {
        fprintf(stderr, "error : invalid bufferSize_s = %g\n", bufferSize_s);
        return false;
    }

    if (bufferSize_s > kMaxBufferSize_s) {
        fprintf(stderr, "error : invalid record size requested - %g s. max allowed is %g s\n", bufferSize_s, kMaxBufferSize_s);
        return false;
    }

    return true;
}
