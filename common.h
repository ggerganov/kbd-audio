/*! \file common.h
 *  \brief Common types and functions
 *  \author Georgi Gerganov
 */

#pragma once

#include <map>
#include <string>
#include <tuple>
#include <vector>

// types

template<typename T>
struct stWaveformView {
    const T * samples = nullptr;
    int64_t n         = 0;
};

template<typename T> using TWaveformT = std::vector<T>;

template<typename T>
stWaveformView<T> getView(const TWaveformT<T> & waveform, int64_t idx) {
    return { waveform.data() + idx, (int64_t) waveform.size() - idx };
}

template<typename T>
stWaveformView<T> getView(const TWaveformT<T> & waveform, int64_t idx, int64_t len) {
    return { waveform.data() + idx, len };
}

using TConfidence = float;
using TValueCC = double;
using TOffset = int64_t;

using TSampleF = float;
using TSampleI16 = int16_t;
using TSampleI32 = int32_t;

using TKey = int;
using TKeyConfidenceMap = std::map<TKey, TConfidence>;

// - i16 samples

using TWaveformI16          = TWaveformT<TSampleI16>;
using TWaveformViewI16      = stWaveformView<TSampleI16>;

// - i32 samples

using TWaveformI32          = TWaveformT<TSampleI32>;
using TWaveformViewI32      = stWaveformView<TSampleI32>;

// - float samples

using TKeyWaveformF = std::vector<TSampleF>;
using TKeyHistoryF = std::vector<TKeyWaveformF>;

// helpers

std::map<std::string, std::string> parseCmdArguments(int argc, char ** argv);

std::tuple<double, double> calcSum(const TKeyWaveformF & waveform, int is0, int is1);

TValueCC calcCC(
    const TKeyWaveformF & waveform0,
    const TKeyWaveformF & waveform1,
    double sum0, double sum02,
    int is00, int is0, int is1);

std::tuple<TValueCC, TOffset> findBestCC(
    const TKeyWaveformF & waveform0,
    const TKeyWaveformF & waveform1,
    int is0, int is1,
    int alignWindow);
