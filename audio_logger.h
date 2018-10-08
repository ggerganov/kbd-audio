/*! \file audio_logger.h
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#pragma once

#include <memory>
#include <array>
#include <vector>
#include <functional>

class AudioLogger {
    public:
        constexpr static int64_t getBufferSize_frames(int64_t sampleRate, float bufferSize_s) {
            return (bufferSize_s*sampleRate)/kSamplesPerFrame;
        }

        static constexpr int64_t kSamplesPerFrame = 512;
        static constexpr int64_t kMaxSampleRate = 96000;
        static constexpr float    kMaxBufferSize_s = 1.000f;

        using Sample = float;
        using Frame = std::array<Sample, kSamplesPerFrame>;
        using Record = std::vector<Frame>;
        using Callback = std::function<void(const Record & frames)>;

        AudioLogger();
        ~AudioLogger();

        bool install(int64_t sampleRate, Callback callback);
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
