/*! \file audio_logger.h
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#pragma once

#include <memory>
#include <array>
#include <list>
#include <functional>

class AudioLogger {
    public:
        static constexpr uint64_t kSampleRate = 44100;
        static constexpr uint64_t kSamplesPerFrame = 1024;
        static constexpr float kBufferSize_s = 0.100;
        static constexpr uint64_t kBufferSize_frames = kBufferSize_s*kSampleRate/kSamplesPerFrame;

        using Sample = float;
        using Frame = std::array<Sample, kSamplesPerFrame>;
        using Record = std::list<Frame>;
        using Callback = std::function<void(const Record & frames)>;

        AudioLogger();
        ~AudioLogger();

        bool install(Callback callback);
        bool addFrame(const Sample * stream);
        bool record();

    private:
        struct Data;
        std::unique_ptr<Data> data_;
        Data & getData() { return *data_; }
};
