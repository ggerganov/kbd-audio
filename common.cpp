/*! \file common.cpp
 *  \brief Enter description here.
 */

#include "common.h"
#include "constants.h"

#include <cstring>
#include <cmath>
#include <thread>
#include <mutex>
#include <fstream>
#include <deque>

namespace {
template <typename TSampleInput, typename TSample>
    bool readWaveform(std::ifstream & fin, TWaveformT<TSample> & res, int32_t offset, std::streamsize size) {
        if (std::is_same<TSample, TSampleI16>::value) {
            std::vector<TSampleInput> buf(size/sizeof(TSampleInput));
            res.resize(offset + size/sizeof(TSampleInput));
            fin.read((char *)(buf.data()), size);
            double amax = 0.0f;
            for (auto i = 0; i < (int) buf.size(); ++i) if (std::abs(buf[i]) > amax) amax = std::abs(buf[i]);
            for (auto i = 0; i < (int) buf.size(); ++i) res[offset + i] = std::round(std::numeric_limits<TSampleI16>::max()*(buf[i]/amax));
        } else if (std::is_same<TSample, TSampleF>::value) {
            res.resize(offset + size/sizeof(TSample));
            fin.read((char *)(res.data() + offset), size);
        } else {
            return false;
        }
        return true;
    }
}

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

template <typename TSampleInput, typename TSample>
bool readFromFile(const std::string & fname, TWaveformT<TSample> & res) {
    std::ifstream fin(fname, std::ios::binary | std::ios::ate);
    if (fin.good() == false) {
        return false;
    }

    {
        std::streamsize size = fin.tellg();
        fin.seekg(0, std::ios::beg);

        static_assert(std::is_same<TSampleInput, TSampleF>::value, "TSampleInput not supported");
        static_assert(std::is_same<TSample, TSampleF>::value ||
                      std::is_same<TSample, TSampleI16>::value, "TSample not supported");

        if (readWaveform<TSampleInput>(fin, res, 0, size) == false) {
            return false;
        }
    }

    fin.close();

    return true;
}

template bool readFromFile<TSampleF, TSampleI16>(const std::string & fname, TWaveformT<TSampleI16> & res);

template <typename TSampleInput, typename TSample>
bool readFromFile(const std::string & fname, TWaveformT<TSample> & res, TTrainKeys & trainKeys, int32_t & bufferSize_frames) {
    trainKeys.clear();

    std::ifstream fin(fname, std::ios::binary);
    if (fin.good() == false) {
        return false;
    }

    bufferSize_frames = 1;
    fin.read((char *)(&bufferSize_frames), sizeof(bufferSize_frames));

    {
        static_assert(std::is_same<TSampleInput, TSampleF>::value, "TSampleInput not supported");
        static_assert(std::is_same<TSample, TSampleF>::value ||
                      std::is_same<TSample, TSampleI16>::value, "TSample not supported");

        int32_t offset = 0;
        std::streamsize size = bufferSize_frames*kSamplesPerFrame*sizeof(TSampleInput);
        while (true) {
            TKey keyPressed = 0;
            fin.read((char *)(&keyPressed), sizeof(keyPressed));
            if (fin.eof()) break;
            trainKeys.push_back(keyPressed);

            if (readWaveform<TSampleInput>(fin, res, offset, size) == false) {
                return false;
            }

            offset += size/sizeof(TSampleInput);
            if (fin.eof()) break;
        }
    }

    fin.close();

    return true;
}

template bool readFromFile<TSampleF, TSampleI16>(const std::string & fname, TWaveformT<TSampleI16> & res, TTrainKeys & trainKeys, int32_t & bufferSize_frames);

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

    auto samples = waveform.samples;
    auto n       = waveform.n;

    for (int is = 0; is < n; ++is) {
        auto a0 = samples[is];
        sum += a0;
        sum2 += a0*a0;
    }

    return std::tuple<int64_t, int64_t>(sum, sum2);
}

// calcSum : specializations

template std::tuple<int64_t, int64_t> calcSum<TSampleI16>(const TWaveformViewT<TSampleI16> & waveform);

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

    auto samples0 = waveform0.samples;
    auto n0       = waveform0.n;

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

template TValueCC calcCC<TSampleI16>(
    const TWaveformViewT<TSampleI16> & waveform0,
    const TWaveformViewT<TSampleI16> & waveform1,
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

    //auto samples0 = waveform0.samples;
    auto n0       = waveform0.n;

    auto samples1 = waveform1.samples;

#ifdef MY_DEBUG
    auto n1 = waveform1.n;
    if (n0 + 2*alignWindow != n1) {
        printf("BUG 924830jm92, n0 = %d, a = %d\n", (int) n0, (int) alignWindow);
    }
#endif

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

template std::tuple<TValueCC, TOffset> findBestCC<TSampleI16>(
    const TWaveformViewT<TSampleI16> & waveform0,
    const TWaveformViewT<TSampleI16> & waveform1,
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
        TSimilarityMap & res) {
    int nPresses = keyPresses.size();

    int w = keyPressWidth_samples;
    int a = alignWindow_samples;

    res.clear();
    res.resize(nPresses);
    for (auto & x : res) x.resize(nPresses);

    for (int i = 0; i < nPresses; ++i) {
        res[i][i].cc = 1.0f;
        res[i][i].offset = 0;

        auto & waveform0 = keyPresses[i].waveform;
        auto & pos0      = keyPresses[i].pos;
        auto & avgcc     = keyPresses[i].ccAvg;

        auto samples0 = waveform0.samples;
        //auto n0       = waveform0.n;

        for (int j = 0; j < nPresses; ++j) {
            if (i == j) continue;

            auto waveform1 = keyPresses[j].waveform;
            auto pos1      = keyPresses[j].pos;

            auto samples1 = waveform1.samples;
            auto ret = findBestCC(TWaveformViewT<T> { samples0 + pos0 + offsetFromPeak_samples,     2*w },
                                  TWaveformViewT<T> { samples1 + pos1 + offsetFromPeak_samples - a, 2*w + 2*a }, a);
            auto bestcc     = std::get<0>(ret);
            auto bestoffset = std::get<1>(ret);

            res[i][j].cc = bestcc;
            res[i][j].offset = bestoffset;

            avgcc += bestcc;
        }
        avgcc /= (nPresses - 1);
    }

    return true;
}

template bool calculateSimilartyMap<TSampleI16>(
        const int32_t keyPressWidth_samples,
        const int32_t alignWindow_samples,
        const int32_t offsetFromPeak_samples,
        TKeyPressCollectionT<TSampleI16> & keyPresses,
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
        int historySize) {
    res.clear();
    waveformThreshold.resize(waveform.n);

    int rbBegin = 0;
    double rbAverage = 0.0;
    std::vector<double> rbSamples(8*historySize, 0.0);

    int k = historySize;
    std::deque<int64_t> que(k);

    auto samples = waveform.samples;
    auto n       = waveform.n;

    TWaveformT<T> waveformAbs(n);
    for (int64_t i = 0; i < n; ++i) {
        waveformAbs[i] = std::abs(samples[i]);
    }

    for (int64_t i = 0; i < n; ++i) {
        {
            int64_t ii = i - k/2;
            if (ii >= 0) {
                rbAverage *= rbSamples.size();
                rbAverage -= rbSamples[rbBegin];
                double acur = waveformAbs[i];
                rbSamples[rbBegin] = acur;
                rbAverage += acur;
                rbAverage /= rbSamples.size();
                if (++rbBegin >= (int) rbSamples.size()) {
                    rbBegin = 0;
                }
            }
        }

        if (i < k) {
            while((!que.empty()) && waveformAbs[i] >= waveformAbs[que.back()]) {
                que.pop_back();
            }
            que.push_back(i);
        } else {
            while((!que.empty()) && que.front() <= i - k) {
                que.pop_front();
            }

            while((!que.empty()) && waveformAbs[i] >= waveformAbs[que.back()]) {
                que.pop_back();
            }

            que.push_back(i);

            int64_t itest = i - k/2;
            if (itest >= 2*k && itest < n - 2*k && que.front() == itest) {
                double acur = waveformAbs[itest];
                if (acur > thresholdBackground*rbAverage){
                    res.emplace_back(TKeyPressDataT<T> { std::move(waveform), itest, 0.0, -1, -1, '?' });
                }
            }
            waveformThreshold[itest] = waveformAbs[que.front()];
        }
    }

    return true;
}

template bool findKeyPresses<TSampleI16>(
        const TWaveformViewT<TSampleI16> & waveform,
        TKeyPressCollectionT<TSampleI16> & res,
        TWaveformT<TSampleI16> & waveformThreshold,
        double thresholdBackground,
        int historySize);
