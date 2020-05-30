/*! \file common.h
 *  \brief Common types and functions
 *  \author Georgi Gerganov
 */

#pragma once

#include <map>
#include <string>
#include <cstring>
#include <tuple>
#include <cmath>
#include <thread>
#include <mutex>
#include <vector>

// types

using TConfidence = float;
using TValueCC = double;
using TOffset = int64_t;

using TSampleF = float;
using TSampleI16 = int16_t;
using TSampleI32 = int32_t;

using TKey = int;
using TKeyConfidenceMap = std::map<TKey, TConfidence>;

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
