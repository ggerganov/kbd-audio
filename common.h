/*! \file common.h
 *  \brief Common types and functions
 *  \author Georgi Gerganov
 */

#pragma once

#include "audio_logger.h"

#include <map>
#include <string>
#include <cstring>
#include <tuple>
#include <cmath>
#include <thread>
#include <mutex>

// types

using TSum = double;
using TSum2 = double;
using TConfidence = float;
using TValueCC = double;
using TOffset = int;

using TKey = int;
using TKeyWaveform = std::vector<AudioLogger::Sample>;
using TKeyHistory = std::vector<TKeyWaveform>;
using TKeyConfidenceMap = std::map<TKey, TConfidence>;

// helpers

static std::map<std::string, std::string> parseCmdArguments(int argc, char ** argv) {
    int last = argc;
    std::map<std::string, std::string> res;
    for (int i = 1; i < last; ++i) {
        if (argv[i][0] == '-') {
            if (strlen(argv[i]) > 1) {
                res[std::string(1, argv[i][1])] = strlen(argv[i]) > 2 ? argv[i] + 2 : "";
            }
        }
    }

    return res;
}

static std::tuple<TSum, TSum2> calcSum(const TKeyWaveform & waveform, int is0, int is1) {
    TSum sum = 0.0f;
    TSum2 sum2 = 0.0f;
    for (int is = is0; is < is1; ++is) {
        auto a0 = waveform[is];
        sum += a0;
        sum2 += a0*a0;
    }

    return std::tuple<TSum, TSum2>(sum, sum2);
}

static TValueCC calcCC(
    const TKeyWaveform & waveform0,
    const TKeyWaveform & waveform1,
    TSum sum0, TSum2 sum02,
    int is00, int is0, int is1) {
    TValueCC cc = -1.0f;

    TSum sum1 = 0.0f;
    TSum2 sum12 = 0.0f;
    TSum2 sum01 = 0.0f;
    for (int is = 0; is < is1 - is0; ++is) {
        auto a0 = waveform0[is00 + is];
        auto a1 = waveform1[is0 + is];

#ifdef MY_DEBUG
        if (is00 + is < 0 || is00 + is >= waveform0.size()) printf("BUG 0\n");
        if (is0 + is < 0 || is0 + is >= waveform1.size()) {
            printf("BUG 1\n");
            printf("%d %d %d\n", is0, is, (int) waveform1.size());
        }
#endif

        sum1 += a1;
        sum12 += a1*a1;
        sum01 += a0*a1;
    }

    int ncc = (is1 - is0);
    {
        double nom = sum01*ncc - sum0*sum1;
        double den2a = sum02*ncc - sum0*sum0;
        double den2b = sum12*ncc - sum1*sum1;
        cc = (nom)/(sqrt(den2a*den2b));
    }

    return cc;
}

std::tuple<TValueCC, TOffset> findBestCC(
    const TKeyWaveform & waveform0,
    const TKeyWaveform & waveform1,
    int is0, int is1,
    int alignWindow) {
    TOffset besto = -1;
    TValueCC bestcc = -1.0f;

    int is00 = waveform0.size()/2 - (is1 - is0)/2;
    //auto [sum0, sum02] = calcSum(waveform0, is00, is00 + is1 - is0);
    auto ret = calcSum(waveform0, is00, is00 + is1 - is0);
    auto sum0  = std::get<0>(ret);
    auto sum02 = std::get<1>(ret);

#ifdef __EMSCRIPTEN__
    TOffset cbesto = -1;
    TValueCC cbestcc = -1.0f;

    for (int o = -alignWindow; o < alignWindow; ++o) {
        auto cc = calcCC(waveform0, waveform1, sum0, sum02, is00, is0 + o, is1 + o);
        if (cc > cbestcc) {
            cbesto = o;
            cbestcc = cc;
        }
    }

    if (cbestcc > bestcc) {
        bestcc = cbestcc;
        besto = cbesto;
    }
#else
    int nWorkers = std::min(4u, std::thread::hardware_concurrency());
    std::mutex mutex;
    std::vector<std::thread> workers(nWorkers);
    for (int i = 0; i < workers.size(); ++i) {
        auto & worker = workers[i];
        worker = std::thread([&, sum0 = sum0, sum02 = sum02, i]() {
            TOffset cbesto = -1;
            TValueCC cbestcc = -1.0f;

            for (int o = -alignWindow + i; o < alignWindow; o += nWorkers) {
                auto cc = calcCC(waveform0, waveform1, sum0, sum02, is00, is0 + o, is1 + o);
                if (cc > cbestcc) {
                    cbesto = o;
                    cbestcc = cc;
                }
            }

            {
                std::lock_guard<std::mutex> lock(mutex);
                if (cbestcc > bestcc) {
                    bestcc = cbestcc;
                    besto = cbesto;
                }
            }
        });
    }
    for (auto & worker : workers) worker.join();
#endif

    return std::tuple<TValueCC, TOffset>(bestcc, besto);
}
