/*! \file audio_logger.h
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#pragma once

#include "constants.h"

#include <memory>
#include <array>
#include <vector>
#include <functional>

class AudioLogger {
    public:
        using Sample = float;
        using Frame = std::array<Sample, kSamplesPerFrame>;
        using Record = std::vector<Frame>;
        using Callback = std::function<void(const Record & frames)>;

        struct Parameters {
            int64_t sampleRate;
            Callback callback;
            int captureId = 0;
            int nChannels = 0;

            enum EFilter {
                None = 0,
                FirstOrderHighPass,
                SecondOrderButterworthHighPass,
            };

            EFilter filter = SecondOrderButterworthHighPass;

            float freqCutoff = 1000.0;
        };

        AudioLogger();
        ~AudioLogger();

        bool install(Parameters && parameters);
        bool terminate();
        bool addFrame(const Sample * stream);
        bool record(float bufferSize_s);
        bool recordSym(float bufferSize_s);

        bool pause();
        bool resume();

    private:
        struct Data;
        std::unique_ptr<Data> data_;
        Data & getData() { return *data_; }
};
