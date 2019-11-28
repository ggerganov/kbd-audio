/*! \file constants.h
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#pragma once

#include <map>
#include <array>

static constexpr int64_t kSamplesPerFrame = 512;
static constexpr int64_t kMaxSampleRate = 96000;
static constexpr float   kMaxBufferSize_s = 1.000f;

static constexpr int32_t ceil_const(float num) {
    return (static_cast<float>(static_cast<int32_t>(num)) == num)
        ? static_cast<int32_t>(num)
        : static_cast<int32_t>(num) + ((num > 0) ? 1 : 0);
}

// odd number of frames longer than bufferSize_s
static constexpr int32_t getBufferSize_frames(int64_t sampleRate, float bufferSize_s) {
    return ceil_const(float(bufferSize_s*sampleRate)/kSamplesPerFrame)/2 + 1;
}

static constexpr int64_t kSampleRate = 24000;
static constexpr float kTrainBufferSize_s = 0.175f;
static constexpr float kPredictBufferSize_s = 0.200f;
static constexpr int32_t kTrainBufferSize_frames = getBufferSize_frames(kSampleRate, kTrainBufferSize_s);
static constexpr int32_t kPredictBufferSize_frames = getBufferSize_frames(kSampleRate, kPredictBufferSize_s);
static constexpr int64_t kSamplesPerWaveform = kSamplesPerFrame*kTrainBufferSize_frames;

static constexpr uint64_t kBkgrRingBufferSize = 4*1024;
static constexpr int64_t kBkgrStep_samples = 1;
static constexpr int64_t kKeyDuration_samples = 0.005f*kSampleRate;

static const std::array<int32_t, 256> kCharToInt = {
    /* { 0,   */   0 /* } */,
    /* { 1,   */   0 /* } */,
    /* { 2,   */   0 /* } */,
    /* { 3,   */   0 /* } */,
    /* { 4,   */   0 /* } */,
    /* { 5,   */   0 /* } */,
    /* { 6,   */   0 /* } */,
    /* { 7,   */   0 /* } */,
    /* { 8,   */   0 /* } */,
    /* { 9,   */   0 /* } */,
    /* { 10,  */   0 /* } */,
    /* { 11,  */   0 /* } */,
    /* { 12,  */   0 /* } */,
    /* { 13,  */   0 /* } */,
    /* { 14,  */   0 /* } */,
    /* { 15,  */   0 /* } */,
    /* { 16,  */   0 /* } */,
    /* { 17,  */   0 /* } */,
    /* { 18,  */   0 /* } */,
    /* { 19,  */   0 /* } */,
    /* { 20,  */   0 /* } */,
    /* { 21,  */   0 /* } */,
    /* { 22,  */   0 /* } */,
    /* { 23,  */   0 /* } */,
    /* { 24,  */   0 /* } */,
    /* { 25,  */   0 /* } */,
    /* { 26,  */   0 /* } */,
    /* { 27,  */   0 /* } */,
    /* { 28,  */   0 /* } */,
    /* { 29,  */   0 /* } */,
    /* { 30,  */   0 /* } */,
    /* { 31,  */   0 /* } */,
    /* { 32,  */   0 /* } */,
    /* { 33,  */   0 /* } */,
    /* { 34,  */   0 /* } */,
    /* { 35,  */   0 /* } */,
    /* { 36,  */   0 /* } */,
    /* { 37,  */   0 /* } */,
    /* { 38,  */   0 /* } */,
    /* { 39,  */   0 /* } */,
    /* { 40,  */   0 /* } */,
    /* { 41,  */   0 /* } */,
    /* { 42,  */   0 /* } */,
    /* { 43,  */   0 /* } */,
    /* { 44,  */   0 /* } */,
    /* { 45,  */   0 /* } */,
    /* { 46,  */   0 /* } */,
    /* { 47,  */   0 /* } */,
    /* { 48,  */   0 /* } */,
    /* { 49,  */   0 /* } */,
    /* { 50,  */   0 /* } */,
    /* { 51,  */   0 /* } */,
    /* { 52,  */   0 /* } */,
    /* { 53,  */   0 /* } */,
    /* { 54,  */   0 /* } */,
    /* { 55,  */   0 /* } */,
    /* { 56,  */   0 /* } */,
    /* { 57,  */   0 /* } */,
    /* { 58,  */   0 /* } */,
    /* { 59,  */   0 /* } */,
    /* { 60,  */   0 /* } */,
    /* { 61,  */   0 /* } */,
    /* { 62,  */   0 /* } */,
    /* { 63,  */   0 /* } */,
    /* { 64,  */   0 /* } */,
    /* { 65,  */   1 /* } */,
    /* { 66,  */   2 /* } */,
    /* { 67,  */   3 /* } */,
    /* { 68,  */   4 /* } */,
    /* { 69,  */   5 /* } */,
    /* { 70,  */   6 /* } */,
    /* { 71,  */   7 /* } */,
    /* { 72,  */   8 /* } */,
    /* { 73,  */   9 /* } */,
    /* { 74,  */  10 /* } */,
    /* { 75,  */  11 /* } */,
    /* { 76,  */  12 /* } */,
    /* { 77,  */  13 /* } */,
    /* { 78,  */  14 /* } */,
    /* { 79,  */  15 /* } */,
    /* { 80,  */  16 /* } */,
    /* { 81,  */  17 /* } */,
    /* { 82,  */  18 /* } */,
    /* { 83,  */  19 /* } */,
    /* { 84,  */  20 /* } */,
    /* { 85,  */  21 /* } */,
    /* { 86,  */  22 /* } */,
    /* { 87,  */  23 /* } */,
    /* { 88,  */  24 /* } */,
    /* { 89,  */  25 /* } */,
    /* { 90,  */  26 /* } */,
    /* { 91,  */   0 /* } */,
    /* { 92,  */   0 /* } */,
    /* { 93,  */   0 /* } */,
    /* { 94,  */   0 /* } */,
    /* { 95,  */   0 /* } */,
    /* { 96,  */   0 /* } */,
    /* { 97,  */   1 /* } */,
    /* { 98,  */   2 /* } */,
    /* { 99,  */   3 /* } */,
    /* { 100, */   4 /* } */,
    /* { 101, */   5 /* } */,
    /* { 102, */   6 /* } */,
    /* { 103, */   7 /* } */,
    /* { 104, */   8 /* } */,
    /* { 105, */   9 /* } */,
    /* { 106, */  10 /* } */,
    /* { 107, */  11 /* } */,
    /* { 108, */  12 /* } */,
    /* { 109, */  13 /* } */,
    /* { 110, */  14 /* } */,
    /* { 111, */  15 /* } */,
    /* { 112, */  16 /* } */,
    /* { 113, */  17 /* } */,
    /* { 114, */  18 /* } */,
    /* { 115, */  19 /* } */,
    /* { 116, */  20 /* } */,
    /* { 117, */  21 /* } */,
    /* { 118, */  22 /* } */,
    /* { 119, */  23 /* } */,
    /* { 120, */  24 /* } */,
    /* { 121, */  25 /* } */,
    /* { 122, */  26 /* } */,
    /* { 123, */   0 /* } */,
    /* { 124, */   0 /* } */,
    /* { 125, */   0 /* } */,
    /* { 126, */   0 /* } */,
    /* { 127, */   0 /* } */,
    /* { 128, */   0 /* } */,
    /* { 129, */   0 /* } */,
    /* { 130, */   0 /* } */,
    /* { 131, */   0 /* } */,
    /* { 132, */   0 /* } */,
    /* { 133, */   0 /* } */,
    /* { 134, */   0 /* } */,
    /* { 135, */   0 /* } */,
    /* { 136, */   0 /* } */,
    /* { 137, */   0 /* } */,
    /* { 138, */   0 /* } */,
    /* { 139, */   0 /* } */,
    /* { 140, */   0 /* } */,
    /* { 141, */   0 /* } */,
    /* { 142, */   0 /* } */,
    /* { 143, */   0 /* } */,
    /* { 144, */   0 /* } */,
    /* { 145, */   0 /* } */,
    /* { 146, */   0 /* } */,
    /* { 147, */   0 /* } */,
    /* { 148, */   0 /* } */,
    /* { 149, */   0 /* } */,
    /* { 150, */   0 /* } */,
    /* { 151, */   0 /* } */,
    /* { 152, */   0 /* } */,
    /* { 153, */   0 /* } */,
    /* { 154, */   0 /* } */,
    /* { 155, */   0 /* } */,
    /* { 156, */   0 /* } */,
    /* { 157, */   0 /* } */,
    /* { 158, */   0 /* } */,
    /* { 159, */   0 /* } */,
    /* { 160, */   0 /* } */,
    /* { 161, */   0 /* } */,
    /* { 162, */   0 /* } */,
    /* { 163, */   0 /* } */,
    /* { 164, */   0 /* } */,
    /* { 165, */   0 /* } */,
    /* { 166, */   0 /* } */,
    /* { 167, */   0 /* } */,
    /* { 168, */   0 /* } */,
    /* { 169, */   0 /* } */,
    /* { 170, */   0 /* } */,
    /* { 171, */   0 /* } */,
    /* { 172, */   0 /* } */,
    /* { 173, */   0 /* } */,
    /* { 174, */   0 /* } */,
    /* { 175, */   0 /* } */,
    /* { 176, */   0 /* } */,
    /* { 177, */   0 /* } */,
    /* { 178, */   0 /* } */,
    /* { 179, */   0 /* } */,
    /* { 180, */   0 /* } */,
    /* { 181, */   0 /* } */,
    /* { 182, */   0 /* } */,
    /* { 183, */   0 /* } */,
    /* { 184, */   0 /* } */,
    /* { 185, */   0 /* } */,
    /* { 186, */   0 /* } */,
    /* { 187, */   0 /* } */,
    /* { 188, */   0 /* } */,
    /* { 189, */   0 /* } */,
    /* { 190, */   0 /* } */,
    /* { 191, */   0 /* } */,
    /* { 192, */   0 /* } */,
    /* { 193, */   0 /* } */,
    /* { 194, */   0 /* } */,
    /* { 195, */   0 /* } */,
    /* { 196, */   0 /* } */,
    /* { 197, */   0 /* } */,
    /* { 198, */   0 /* } */,
    /* { 199, */   0 /* } */,
    /* { 200, */   0 /* } */,
    /* { 201, */   0 /* } */,
    /* { 202, */   0 /* } */,
    /* { 203, */   0 /* } */,
    /* { 204, */   0 /* } */,
    /* { 205, */   0 /* } */,
    /* { 206, */   0 /* } */,
    /* { 207, */   0 /* } */,
    /* { 208, */   0 /* } */,
    /* { 209, */   0 /* } */,
    /* { 210, */   0 /* } */,
    /* { 211, */   0 /* } */,
    /* { 212, */   0 /* } */,
    /* { 213, */   0 /* } */,
    /* { 214, */   0 /* } */,
    /* { 215, */   0 /* } */,
    /* { 216, */   0 /* } */,
    /* { 217, */   0 /* } */,
    /* { 218, */   0 /* } */,
    /* { 219, */   0 /* } */,
    /* { 220, */   0 /* } */,
    /* { 221, */   0 /* } */,
    /* { 222, */   0 /* } */,
    /* { 223, */   0 /* } */,
    /* { 224, */   0 /* } */,
    /* { 225, */   0 /* } */,
    /* { 226, */   0 /* } */,
    /* { 227, */   0 /* } */,
    /* { 228, */   0 /* } */,
    /* { 229, */   0 /* } */,
    /* { 230, */   0 /* } */,
    /* { 231, */   0 /* } */,
    /* { 232, */   0 /* } */,
    /* { 233, */   0 /* } */,
    /* { 234, */   0 /* } */,
    /* { 235, */   0 /* } */,
    /* { 236, */   0 /* } */,
    /* { 237, */   0 /* } */,
    /* { 238, */   0 /* } */,
    /* { 239, */   0 /* } */,
    /* { 240, */   0 /* } */,
    /* { 241, */   0 /* } */,
    /* { 242, */   0 /* } */,
    /* { 243, */   0 /* } */,
    /* { 244, */   0 /* } */,
    /* { 245, */   0 /* } */,
    /* { 246, */   0 /* } */,
    /* { 247, */   0 /* } */,
    /* { 248, */   0 /* } */,
    /* { 249, */   0 /* } */,
    /* { 250, */   0 /* } */,
    /* { 251, */   0 /* } */,
    /* { 252, */   0 /* } */,
    /* { 253, */   0 /* } */,
    /* { 254, */   0 /* } */,
    /* { 255, */   0 /* } */,
};

static const std::map<int, const char *> kKeyText = {
    { -1,  "?" },
    { 0,   "?" },
    { 1,   "?" },
    { 2,   "?" },
    { 3,   "?" },
    { 4,   "?" },
    { 5,   "?" },
    { 6,   "?" },
    { 7,   "?" },
    { 8,   "?" },
    { 9,   "?" },
    { 10,  "[enter]" },
    { 11,  "?" },
    { 12,  "?" },
    { 13,  "?" },
    { 14,  "?" },
    { 15,  "?" },
    { 16,  "?" },
    { 17,  "?" },
    { 18,  "?" },
    { 19,  "?" },
    { 20,  "?" },
    { 21,  "?" },
    { 22,  "?" },
    { 23,  "?" },
    { 24,  "?" },
    { 25,  "?" },
    { 26,  "?" },
    { 27,  "?" },
    { 28,  "?" },
    { 29,  "?" },
    { 30,  "?" },
    { 31,  "?" },
    { 32,  "[space]" },
    { 33,  "?" },
    { 34,  "?" },
    { 35,  "?" },
    { 36,  "?" },
    { 37,  "?" },
    { 38,  "?" },
    { 39,  "'" },
    { 40,  "?" },
    { 41,  "?" },
    { 42,  "?" },
    { 43,  "?" },
    { 44,  "," },
    { 45,  "-" },
    { 46,  "." },
    { 47,  "/" },
    { 48,  "0" },
    { 49,  "1" },
    { 50,  "2" },
    { 51,  "3" },
    { 52,  "4" },
    { 53,  "5" },
    { 54,  "6" },
    { 55,  "7" },
    { 56,  "8" },
    { 57,  "9" },
    { 58,  "?" },
    { 59,  ";" },
    { 60,  "?" },
    { 61,  "=" },
    { 62,  "?" },
    { 63,  "?" },
    { 64,  "?" },
    { 65,  "?" },
    { 66,  "?" },
    { 67,  "?" },
    { 68,  "?" },
    { 69,  "?" },
    { 70,  "?" },
    { 71,  "?" },
    { 72,  "?" },
    { 73,  "?" },
    { 74,  "?" },
    { 75,  "?" },
    { 76,  "?" },
    { 77,  "?" },
    { 78,  "?" },
    { 79,  "?" },
    { 80,  "?" },
    { 81,  "?" },
    { 82,  "?" },
    { 83,  "?" },
    { 84,  "?" },
    { 85,  "?" },
    { 86,  "?" },
    { 87,  "?" },
    { 88,  "?" },
    { 89,  "?" },
    { 90,  "?" },
    { 91,  "[" },
    { 92,  "\\" },
    { 93,  "]" },
    { 94,  "?" },
    { 95,  "?" },
    { 96,  "`" },
    { 97,  "a" },
    { 98,  "b" },
    { 99,  "c" },
    { 100, "d" },
    { 101, "e" },
    { 102, "f" },
    { 103, "g" },
    { 104, "h" },
    { 105, "i" },
    { 106, "j" },
    { 107, "k" },
    { 108, "l" },
    { 109, "m" },
    { 110, "n" },
    { 111, "o" },
    { 112, "p" },
    { 113, "q" },
    { 114, "r" },
    { 115, "s" },
    { 116, "t" },
    { 117, "u" },
    { 118, "v" },
    { 119, "w" },
    { 120, "x" },
    { 121, "y" },
    { 122, "z" },
    { 123, "?" },
    { 124, "?" },
    { 125, "?" },
    { 126, "?" },
    { 127, "[<-]" },
    { 128, "?" },
    { 129, "?" },
    { 130, "?" },
    { 131, "?" },
    { 132, "?" },
    { 133, "?" },
    { 134, "?" },
    { 135, "?" },
    { 136, "?" },
    { 137, "?" },
    { 138, "?" },
    { 139, "?" },
    { 140, "?" },
    { 141, "?" },
    { 142, "?" },
    { 143, "?" },
    { 144, "?" },
    { 145, "?" },
    { 146, "?" },
    { 147, "?" },
    { 148, "?" },
    { 149, "?" },
    { 150, "?" },
    { 151, "?" },
    { 152, "?" },
    { 153, "?" },
    { 154, "?" },
    { 155, "?" },
    { 156, "?" },
    { 157, "?" },
    { 158, "?" },
    { 159, "?" },
    { 160, "?" },
    { 161, "?" },
    { 162, "?" },
    { 163, "?" },
    { 164, "?" },
    { 165, "?" },
    { 166, "?" },
    { 167, "?" },
    { 168, "?" },
    { 169, "?" },
    { 170, "?" },
    { 171, "?" },
    { 172, "?" },
    { 173, "?" },
    { 174, "?" },
    { 175, "?" },
    { 176, "?" },
    { 177, "?" },
    { 178, "?" },
    { 179, "?" },
    { 180, "?" },
    { 181, "?" },
    { 182, "?" },
    { 183, "?" },
    { 184, "?" },
    { 185, "?" },
    { 186, "?" },
    { 187, "?" },
    { 188, "?" },
    { 189, "?" },
    { 190, "?" },
    { 191, "?" },
    { 192, "?" },
    { 193, "?" },
    { 194, "?" },
    { 195, "?" },
    { 196, "?" },
    { 197, "?" },
    { 198, "?" },
    { 199, "?" },
    { 200, "?" },
    { 201, "?" },
    { 202, "?" },
    { 203, "?" },
    { 204, "?" },
    { 205, "?" },
    { 206, "?" },
    { 207, "?" },
    { 208, "?" },
    { 209, "?" },
    { 210, "?" },
    { 211, "?" },
    { 212, "?" },
    { 213, "?" },
    { 214, "?" },
    { 215, "?" },
    { 216, "?" },
    { 217, "?" },
    { 218, "?" },
    { 219, "?" },
    { 220, "?" },
    { 221, "?" },
    { 222, "?" },
    { 223, "?" },
    { 224, "?" },
    { 225, "?" },
    { 226, "?" },
    { 227, "?" },
    { 228, "?" },
    { 229, "?" },
    { 230, "?" },
    { 231, "?" },
    { 232, "?" },
    { 233, "?" },
    { 234, "?" },
    { 235, "?" },
    { 236, "?" },
    { 237, "?" },
    { 238, "?" },
    { 239, "?" },
    { 240, "?" },
    { 241, "?" },
    { 242, "?" },
    { 243, "?" },
    { 244, "?" },
    { 245, "?" },
    { 246, "?" },
    { 247, "?" },
    { 248, "?" },
    { 249, "?" },
    { 250, "?" },
    { 251, "?" },
    { 252, "?" },
    { 253, "?" },
    { 254, "?" },
    { 255, "?" },
};
