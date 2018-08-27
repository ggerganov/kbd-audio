/*! \file audio_logger.h
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#pragma once

#include <memory>

class AudioLogger {
    public:
        static constexpr uint64_t kSamplesPerFrame = 1024;
        static constexpr uint64_t kFramesPerPatch = 1024;

        using Sample = int16_t;
        using Frame = std::array<Sample, kSamplesPerFrame>;
        using Patch = std::array<Frame, kFramesPerPatch>;

        AudioLogger();
        ~AudioLogger();

        bool install();

    private:
        struct Data;
        std::unique_ptr<Data> data_;
        Data & getData() { return *data_; }
};
