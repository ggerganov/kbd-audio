/*! \file common.cpp
 *  \brief Enter description here.
 */

#include "common.h"

#include <cstring>
#include <cmath>
#include <thread>
#include <mutex>

std::map<std::string, std::string> parseCmdArguments(int argc, char ** argv) {
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

//
// calcCC
//

std::tuple<double, double> calcSum(const TKeyWaveformF & waveform, int is0, int is1) {
    double sum = 0.0f;
    double sum2 = 0.0f;
    for (int is = is0; is < is1; ++is) {
        auto a0 = waveform[is];
        sum += a0;
        sum2 += a0*a0;
    }

    return std::tuple<double, double>(sum, sum2);
}

template<typename T>
std::tuple<int64_t, int64_t> calcSum(const TWaveformViewT<T> & waveform) {
    int64_t sum = 0;
    int64_t sum2 = 0;

    //auto [samples, n] = waveform;
    auto samples = waveform.samples;
    auto n       = waveform.n;

    for (int64_t is = 0; is < n; ++is) {
        auto a0 = samples[is];
        sum += a0;
        sum2 += a0*a0;
    }

    return std::tuple<int64_t, int64_t>(sum, sum2);
}

// calcSum : specializations

template std::tuple<int64_t, int64_t> calcSum<TSampleI32>(const TWaveformViewT<TSampleI32> & waveform);

//
// calcCC
//

TValueCC calcCC(
    const TKeyWaveformF & waveform0,
    const TKeyWaveformF & waveform1,
    double sum0, double sum02,
    int is00, int is0, int is1) {
    TValueCC cc = -1.0f;

    double sum1 = 0.0f;
    double sum12 = 0.0f;
    double sum01 = 0.0f;
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

template<typename T>
TValueCC calcCC(
    const TWaveformViewT<T> & waveform0,
    const TWaveformViewT<T> & waveform1,
    int64_t sum0, int64_t sum02) {
    TValueCC cc = -1.0f;

    int64_t sum1 = 0;
    int64_t sum12 = 0;
    int64_t sum01 = 0;

    //auto [samples0, n0] = waveform0;
    auto samples0 = waveform0.samples;
    auto n0       = waveform0.n;

    //auto [samples1, n1] = waveform1;
    auto samples1 = waveform1.samples;
    auto n1       = waveform1.n;

#ifdef MY_DEBUG
    if (n0 != n1) {
        printf("BUG 234f8273\n");
    }
#endif
    auto n = std::min(n0, n1);

    for (int64_t is = 0; is < n; ++is) {
        auto a0 = samples0[is];
        auto a1 = samples1[is];

        sum1 += a1;
        sum12 += a1*a1;
        sum01 += a0*a1;
    }

    {
        double nom = sum01*n - sum0*sum1;
        double den2a = sum02*n - sum0*sum0;
        double den2b = sum12*n - sum1*sum1;
        cc = (nom)/(sqrt(den2a*den2b));
    }

    return cc;
}

// calcCC : specializations

template TValueCC calcCC<TSampleI32>(
    const TWaveformViewT<TSampleI32> & waveform0,
    const TWaveformViewT<TSampleI32> & waveform1,
    int64_t sum0, int64_t sum02);

//
// findBestCC
//

std::tuple<TValueCC, TOffset> findBestCC(
    const TKeyWaveformF & waveform0,
    const TKeyWaveformF & waveform1,
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
    for (int i = 0; i < (int) workers.size(); ++i) {
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

template<typename T>
std::tuple<TValueCC, TOffset> findBestCC(
    const TWaveformViewT<T> & waveform0,
    const TWaveformViewT<T> & waveform1,
    int64_t alignWindow) {
    TValueCC bestcc = -1.0;
    TOffset besto = -1;

    //auto [samples0, n0] = waveform0;
    //auto samples0 = waveform0.samples;
    auto n0       = waveform0.n;

    //auto [samples1, n1] = waveform1;
    auto samples1 = waveform1.samples;

#ifdef MY_DEBUG
    auto n1 = waveform1.n;
    if (n0 + 2*alignWindow != n1) {
        printf("BUG 924830jm92, n0 = %d, a = %d\n", (int) n0, (int) alignWindow);
    }
#endif

    //auto [sum0, sum02] = calcSum(waveform0);
    auto ret = calcSum(waveform0);
    auto sum0  = std::get<0>(ret);
    auto sum02 = std::get<1>(ret);

    for (int o = 0; o < 2*alignWindow; ++o) {
        auto cc = calcCC(waveform0, { samples1 + o, n0 }, sum0, sum02);
        if (cc > bestcc) {
            besto = o - alignWindow;
            bestcc = cc;
        }
    }

    return std::tuple<TValueCC, TOffset>(bestcc, besto);
}

// findBestCC : specializations

template std::tuple<TValueCC, TOffset> findBestCC<TSampleI32>(
    const TWaveformViewT<TSampleI32> & waveform0,
    const TWaveformViewT<TSampleI32> & waveform1,
    int64_t alignWindow);

