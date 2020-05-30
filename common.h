/*! \file common.h
 *  \brief Common types and functions
 *  \author Georgi Gerganov
 */

#pragma once

#include <map>
#include <string>
#include <tuple>
#include <vector>
#include <chrono>

// types

struct stMatch;
template<typename T> struct stWaveformView;
template<typename T> struct stKeyPressData;
template<typename T> struct stKeyPressCollection;

template<typename T> using TWaveformT           = std::vector<T>;
template<typename T> using TWaveformViewT       = stWaveformView<T>;
template<typename T> using TKeyPressDataT       = stKeyPressData<T>;
template<typename T> using TKeyPressCollectionT = stKeyPressCollection<T>;

using TConfidence   = float;
using TValueCC      = double;
using TOffset       = int64_t;
using TClusterId    = int32_t;

using TSampleF      = float;
using TSampleI16    = int16_t;

using TKey              = int32_t;
using TKeyPressPosition = int64_t;
using TKeyConfidenceMap = std::map<TKey, TConfidence>;
using TTrainKeys        = std::vector<TKey>;

using TMatch         = stMatch;
using TSimilarityMap = std::vector<std::vector<TMatch>>;

// - i16 samples

using TWaveformI16              = TWaveformT<TSampleI16>;
using TWaveformViewI16          = TWaveformViewT<TSampleI16>;
using TKeyPressDataI16          = TKeyPressDataT<TSampleI16>;
using TKeyPressCollectionI16    = TKeyPressCollectionT<TSampleI16>;

// - float samples

using TKeyWaveformF = std::vector<TSampleF>;
using TKeyHistoryF = std::vector<TKeyWaveformF>;

// structs
struct stMatch {
    TValueCC    cc      = 0.0;
    TOffset     offset  = 0;
};

template<typename T>
struct stWaveformView {
    const T * samples = nullptr;
    int64_t n         = 0;
};

template<typename T>
struct stKeyPressData {
    TWaveformViewT<T>   waveform;
    TKeyPressPosition   pos         = 0;
    TValueCC            ccAvg       = 0.0;
    TClusterId          cid         = -1;
    TKey                bind        = -1;
    TKey                predicted   = '?';
};

template<typename T>
struct stKeyPressCollection : public std::vector<TKeyPressDataT<T>> {
    int nClusters = 0;
};

// helpers

inline float frand() { return ((float)rand())/RAND_MAX; }

template<typename T>
stWaveformView<T> getView(const TWaveformT<T> & waveform, int64_t idx) {
    return { waveform.data() + idx, (int64_t) waveform.size() - idx };
}

template<typename T>
stWaveformView<T> getView(const TWaveformT<T> & waveform, int64_t idx, int64_t len) {
    return { waveform.data() + idx, len };
}

std::map<std::string, std::string> parseCmdArguments(int argc, char ** argv);

template <typename T>
float toSeconds(T t0, T t1) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()/1000.0f;
}

template <typename TSampleInput, typename TSample>
bool readFromFile(const std::string & fname, TWaveformT<TSample> & res);

template <typename TSampleInput, typename TSample>
bool readFromFile(const std::string & fname, TWaveformT<TSample> & res, TTrainKeys & trainKeys, int32_t & bufferSize_frames);

//
// calcSum
//

std::tuple<double, double> calcSum(const TKeyWaveformF & waveform, int is0, int is1);

template<typename T>
std::tuple<int64_t, int64_t> calcSum(const TWaveformViewT<T> & waveform);

//
// calcCC
//

TValueCC calcCC(
    const TKeyWaveformF & waveform0,
    const TKeyWaveformF & waveform1,
    double sum0, double sum02,
    int is00, int is0, int is1);

template<typename T>
TValueCC calcCC(
    const TWaveformViewT<T> & waveform0,
    const TWaveformViewT<T> & waveform1,
    int64_t sum0, int64_t sum02);

//
// findBestCC
//

std::tuple<TValueCC, TOffset> findBestCC(
    const TKeyWaveformF & waveform0,
    const TKeyWaveformF & waveform1,
    int is0, int is1,
    int alignWindow);

template<typename T>
std::tuple<TValueCC, TOffset> findBestCC(
    const TWaveformViewT<T> & waveform0,
    const TWaveformViewT<T> & waveform1,
    int64_t alignWindow);

//
// calculateSimilarityMap
//

template<typename T>
bool calculateSimilartyMap(
        const int32_t keyPressWidth_samples,
        const int32_t alignWindow_samples,
        const int32_t offsetFromPeak_samples,
        TKeyPressCollectionT<T> & keyPresses,
        TSimilarityMap & res);

//
// findKeyPresses
//

template<typename T>
bool findKeyPresses(
        const TWaveformViewT<T> & waveform,
        TKeyPressCollectionT<T> & res,
        TWaveformT<T> & waveformThreshold,
        double thresholdBackground,
        int historySize);
