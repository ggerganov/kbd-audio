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

        AudioLogger();
        ~AudioLogger();

        bool install(int64_t sampleRate, Callback callback, int captureId = 0);
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
