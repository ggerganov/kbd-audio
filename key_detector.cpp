/*! \file key_detector.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "audio_logger.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <chrono>
#include <thread>

int main(int , const char ** ) {
    constexpr float kBufferSize_s = 0.150f;
    constexpr uint64_t kSampleRate = 96000;
    constexpr uint64_t kRingBufferSize = 16*1024;

    constexpr int bkgrStep_samples = 7;
    constexpr int keyDuration_samples = 0.150f*kSampleRate;

    // rig buffer
    int rbBegin = 0;
    float rbAverage = 0.0f;
    std::array<float, kRingBufferSize> rbSamples;
    rbSamples.fill(0.0f);

    bool doRecord = true;

    int nFrames = getBufferSize_frames(kSampleRate, kBufferSize_s);
    int nFrames2 = std::max(1, nFrames/2);

    AudioLogger audioLogger;
    AudioLogger::Callback cbAudio = [&](const auto & frames) {
        doRecord = true;
        float amax = 0.0f;
        for (int f = 0; f < (int) frames.size(); ++f) {
            for (int s = 0; s < (int) frames[f].size(); s += bkgrStep_samples) {
                rbAverage *= rbSamples.size();
                rbAverage -= rbSamples[rbBegin];
                auto acur = std::abs(frames[f][s]);
                rbSamples[rbBegin] = acur;
                if (acur > amax) amax = acur;
                rbAverage += acur;
                rbAverage /= rbSamples.size();
                if (++rbBegin >= (int) rbSamples.size()) rbBegin = 0;
            }
        }

        int skip_samples = 0;
        int nFrames = frames.size();
        int nFrames2 = std::max(1, nFrames/2);
        for (int f = nFrames2 - nFrames2/2; f <= nFrames2 + nFrames2/2; ++f) {
            for (int s = 0; s < (int) frames[f].size(); ++s) {
                if (s + skip_samples >= (int) frames[f].size()) {
                    skip_samples -= frames[f].size() - s;
                    s += skip_samples;
                    continue;
                } else {
                    s += skip_samples;
                    skip_samples = 0;
                }
                auto acur = frames[f][s];
                if (acur > 10.0f*rbAverage) {
                    skip_samples = keyDuration_samples;
                    printf("Key press detected\n");
                }
            }
        }

        printf("Average = %10.8f, max = %10.8f\n", rbAverage, amax);
    };

    AudioLogger::Parameters parameters;
    parameters.callback = std::move(cbAudio);
    parameters.captureId = 0;
    parameters.nChannels = 1;
    parameters.sampleRate = kSampleRate;
    parameters.freqCutoff_Hz = kFreqCutoff_Hz;

    if (audioLogger.install(std::move(parameters)) == false) {
        fprintf(stderr, "Failed to install audio logger\n");
        return -1;
    }

    while (true) {
        if (doRecord) {
            doRecord = false;
            audioLogger.record(kBufferSize_s, nFrames2/2);
        } else {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
    }

    return 0;
}
