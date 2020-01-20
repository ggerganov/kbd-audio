/*! \file audio_logger.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "audio_logger.h"

#include <SDL.h>
#include <SDL_audio.h>

#include <mutex>
#include <atomic>
#include <algorithm>

#ifndef pi
#define  pi 3.1415926535897932384626433832795
#endif

#ifndef sqrt2
#define sqrt2 (2.0 * 0.707106781186547524401)
#endif

#ifndef sqrt2over2
#define  sqrt2over2  0.707106781186547524401
#endif

namespace {

    // DSP filter
    // ref : https://github.com/dimtass/DSP-Cpp-filters

    struct TFilterCoefficients {
        float a0;
        float a1;
        float a2;
        float b1;
        float b2;
        float c0;
        float d0;

        float xnz1;
        float xnz2;
        float ynz1;
        float ynz2;
    };

    TFilterCoefficients calculateCoefficientsFirstOrderHighPass(int fc, int fs) {
        TFilterCoefficients res;

        float th = 2.0 * pi * fc / fs;
        float g = cos(th) / (1.0 + sin(th));
        res.a0 = (1.0 + g) / 2.0;
        res.a1 = -((1.0 + g) / 2.0);
        res.a2 = 0.0;
        res.b1 = -g;
        res.b2 = 0.0;

        return res;
    }

    TFilterCoefficients calculateCoefficientsSecondOrderButterworthHighPass(int fc, int fs) {
        TFilterCoefficients res;

        float c = tan(pi*fc / fs);
        res.a0 = 1.0 / (1.0 + sqrt2*c + pow(c, 2.0));
        res.a1 = -2.0 * res.a0;
        res.a2 = res.a0;
        res.b1 = 2.0 * res.a0*(pow(c, 2.0) - 1.0);
        res.b2 = res.a0 * (1.0 - sqrt2*c + pow(c, 2.0));

        return res;
    }

    AudioLogger::Sample filterFirstOrderHighPass(TFilterCoefficients & coefficients, AudioLogger::Sample sample) {
        AudioLogger::Sample xn = sample;
        AudioLogger::Sample yn =
            coefficients.a0*xn + coefficients.a1*coefficients.xnz1 + coefficients.a2*coefficients.xnz2 -
            coefficients.b1*coefficients.ynz1 - coefficients.b2*coefficients.xnz2;

        coefficients.xnz2 = coefficients.xnz1;
        coefficients.xnz1 = xn;
        coefficients.xnz2 = coefficients.ynz1;
        coefficients.ynz1 = yn;

        return yn;
    }

    AudioLogger::Sample filterSecondOrderButterworthHighPass(TFilterCoefficients & coefficients, AudioLogger::Sample sample) {
        AudioLogger::Sample xn = sample;
        AudioLogger::Sample yn =
            coefficients.a0*xn + coefficients.a1*coefficients.xnz1 + coefficients.a2*coefficients.xnz2 -
            coefficients.b1*coefficients.ynz1 - coefficients.b2*coefficients.xnz2;

        coefficients.xnz2 = coefficients.xnz1;
        coefficients.xnz1 = xn;
        coefficients.xnz2 = coefficients.ynz1;
        coefficients.ynz1 = yn;

        return yn;
    }

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

    Callback callback = nullptr;

    int64_t sampleRate = kMaxSampleRate;

    int32_t nChannels = 0;
    int32_t sampleSize_bytes = 4;

    int32_t bufferId = 0;
    std::array<Frame, getBufferSize_frames(kMaxSampleRate, kMaxBufferSize_s)> buffer;

    int32_t nRecords = 0;
    std::array<int32_t, kMaxRecords>  nFramesToRecord;
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

    if (SDL_Init(SDL_INIT_AUDIO) < 0) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't initialize SDL: %s\n", SDL_GetError());
        return false;
    }

    int nDevices = SDL_GetNumAudioDevices(SDL_TRUE);
    printf("Found %d capture devices:\n", nDevices);
    for (int i = 0; i < nDevices; i++) {
        printf("    - Capture device #%d: '%s'\n", i, SDL_GetAudioDeviceName(i, SDL_TRUE));
    }

    if (parameters.captureId < 0 || parameters.captureId >= nDevices) {
        printf("Invalid capture device id selected - %d\n", parameters.captureId);
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

    if (obtainedSpec.channels != parameters.nChannels && parameters.nChannels != 0) {
        printf("Warning: obtained number of channels (%d) does not match requested (%d)\n",
               obtainedSpec.channels, parameters.nChannels);
    }

    SDL_PauseAudioDevice(data.deviceIdIn, 0);

    data.sampleRate = parameters.sampleRate;
    data.callback = std::move(parameters.callback);
    data.nChannels = obtainedSpec.channels;
    data.isReady = true;

    switch (parameters.FirstOrderHighPass) {
        case AudioLogger::Parameters::EFilter::None:
            {
            }
            break;
        case AudioLogger::Parameters::EFilter::FirstOrderHighPass:
            {
                data.filterCoefficients = ::calculateCoefficientsFirstOrderHighPass(parameters.freqCutoff, data.sampleRate);
            }
            break;
        case AudioLogger::Parameters::EFilter::SecondOrderButterworthHighPass:
            {
                data.filterCoefficients = ::calculateCoefficientsSecondOrderButterworthHighPass(parameters.freqCutoff, data.sampleRate);
            }
            break;
    };

    data.parameters = parameters;

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

	if ((int) SDL_GetQueuedAudioSize(data.deviceIdIn) > 32*sizeof(float)*kSamplesPerFrame) {
		printf("Queue size: %d\n", SDL_GetQueuedAudioSize(data.deviceIdIn));
		SDL_ClearQueuedAudio(data.deviceIdIn);
	}

    const float norm = 1.0/data.nChannels;

    auto & curFrame = data.buffer[data.bufferId];

    for (int i = 0; i < kSamplesPerFrame; ++i) {
        Sample x = 0;
        for (int j = 0; j < data.nChannels; ++j) {
            x += stream[i*data.nChannels + j];
        }
        curFrame[i] = x*norm;
    }

    switch (data.parameters.FirstOrderHighPass) {
        case AudioLogger::Parameters::EFilter::None:
            {
            }
            break;
        case AudioLogger::Parameters::EFilter::FirstOrderHighPass:
            {
                for (auto & s : curFrame) {
                    s = ::filterFirstOrderHighPass(data.filterCoefficients, s);
                }
            }
            break;
        case AudioLogger::Parameters::EFilter::SecondOrderButterworthHighPass:
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
                if (data.callback) data.callback(record);
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

bool AudioLogger::record(float bufferSize_s) {
    auto & data = getData();

    if (bufferSize_s > kMaxBufferSize_s) {
        fprintf(stderr, "warning :invalid record size requested - %g s. max allowedis %g s\n", bufferSize_s, kMaxBufferSize_s);
        return false;
    }

    auto bufferSize_frames = getBufferSize_frames(data.sampleRate, bufferSize_s);

    std::lock_guard<std::mutex> lock(data.mutex);

    if (data.nRecords == kMaxRecords) {
        fprintf(stderr, "warning : max number of simultaneous records %d reached\n", kMaxRecords);
        return false;
    }

    auto & record = data.records[data.nRecords];

    int nPrev = 3;
    if (record.size() == 0) {
        int fStart = data.bufferId - nPrev;
        if (fStart < 0) fStart += data.buffer.size();
        for (size_t i = 0; i < nPrev; ++i) {
            record.push_back(data.buffer[(fStart + i)%data.buffer.size()]);
        }
    } else {
        fprintf(stderr, "warning : new record requested before last has been processed. should never happen\n");
    }

    data.nFramesToRecord[data.nRecords] = bufferSize_frames - nPrev;

    ++data.nRecords;

    return true;
}

bool AudioLogger::recordSym(float bufferSize_s) {
    auto & data = getData();

    if (bufferSize_s > kMaxBufferSize_s) {
        fprintf(stderr, "warning :invalid record size requested - %g s. max allowedis %g s\n", bufferSize_s, kMaxBufferSize_s);
        return false;
    }

    auto bufferSize_frames = getBufferSize_frames(data.sampleRate, bufferSize_s);

    if (data.nRecords == kMaxRecords) {
        fprintf(stderr, "warning : max number of simultaneous records %d reached\n", kMaxRecords);
        return false;
    }

    std::lock_guard<std::mutex> lock(data.mutex);

    auto & record = data.records[data.nRecords];

    if (record.size() == 0) {
        int fStart = data.bufferId - bufferSize_frames/2;
        if (fStart < 0) fStart += data.buffer.size();
        for (size_t i = 0; i < bufferSize_frames/2; ++i) {
            record.push_back(data.buffer[(fStart + i)%data.buffer.size()]);
        }
    } else {
        fprintf(stderr, "warning : new record requested before last has been processed. should never happen\n");
    }

    data.nFramesToRecord[data.nRecords] = bufferSize_frames/2 + 1;

    ++data.nRecords;

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
