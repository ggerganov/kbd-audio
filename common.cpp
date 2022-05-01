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
#include <algorithm>

#ifndef pi
#define  pi 3.1415926535897932384626433832795
#endif

#ifndef sqrt2
#define sqrt2 (2.0 * 0.707106781186547524401)
#endif

#ifndef sqrt2over2
#define  sqrt2over2  0.707106781186547524401
#endif

namespace {
template <typename TSampleInput, typename TSample>
    bool readWaveform(std::ifstream & fin, TWaveformT<TSample> & res, int32_t offset, std::streamsize size) {
        if (std::is_same<TSample, TSampleI16>::value) {
            std::vector<TSampleInput> buf(size/sizeof(TSampleInput));
            res.resize(offset + size/sizeof(TSampleInput));
            fin.read((char *)(buf.data()), size);
            double amax = calcAbsMax(buf);
            double iamax = amax != 0.0 ? 1.0/amax : 1.0;
            for (auto i = 0; i < (int) buf.size(); ++i) res[offset + i] = std::round(std::numeric_limits<TSampleI16>::max()*(buf[i]*iamax));
        } else if (std::is_same<TSample, TSampleF>::value) {
            res.resize(offset + size/sizeof(TSample));
            fin.read((char *)(res.data() + offset), size);
        } else {
            return false;
        }
        return true;
    }
}

constexpr float iRAND_MAX = 1.0f/float(RAND_MAX);
float frand() { return ((float)rand())*iRAND_MAX; }

float frandGaussian(float mu, float sigma) {
	static const float two_pi = 2.0*3.14159265358979323846;

	thread_local float z1;
	thread_local bool generate;
	generate = !generate;

	if (!generate)
	   return z1 * sigma + mu;

	float u1 = frand();
    float u2 = frand();

    float t = sqrt(-2.0f * log(1.0f - u1));

	float z0 = t*cos(two_pi*u2);
	z1 = t*sin(two_pi*u2);

	return z0 * sigma + mu;
}

uint64_t t_ms() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(); // duh ..
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

template <typename TSampleSrc, typename TSampleDst>
bool convert(const TWaveformT<TSampleSrc> & src, TWaveformT<TSampleDst> & dst) {
    static_assert(std::is_same<TSampleSrc, TSampleDst>::value == false, "Required different sample types");

    static_assert(std::is_same<TSampleSrc, TSampleF>::value, "Source sample type not supported");
    static_assert(std::is_same<TSampleDst, TSampleI16>::value, "Destination sample type not supported");

    dst.resize(src.size());

    double amax = calcAbsMax(src);
    double iamax = amax != 0.0 ? 1.0/amax : 1.0;
    for (auto i = 0; i < (int) src.size(); ++i) dst[i] = std::round(std::numeric_limits<TSampleDst>::max()*(src[i]*iamax));

    return true;
}

template bool convert<TSampleF, TSampleI16>(const TWaveformT<TSampleF> & src, TWaveformT<TSampleI16> & dst);

template <typename TSample>
void filter(TWaveformT<TSample> & waveform, EAudioFilter filterId, float freqCutoff_Hz, int64_t sampleRate) {
    switch (filterId) {
        case EAudioFilter::None:
            {
                return;
            }
            break;
        case EAudioFilter::FirstOrderHighPass:
            {
                auto filterCoefficients = ::calculateCoefficientsFirstOrderHighPass(freqCutoff_Hz, sampleRate);
                for (auto & s : waveform) {
                    s = ::filterFirstOrderHighPass(filterCoefficients, s);
                }
                return;
            }
            break;
        case EAudioFilter::SecondOrderButterworthHighPass:
            {
                auto filterCoefficients = ::calculateCoefficientsSecondOrderButterworthHighPass(freqCutoff_Hz, sampleRate);
                for (auto & s : waveform) {
                    s = ::filterSecondOrderButterworthHighPass(filterCoefficients, s);
                }
                return;
            }
            break;
    }

    fprintf(stderr, "Unknown filter type: %d\n", filterId);
}

template void filter<TSampleF>(TWaveformT<TSampleF> & waveform, EAudioFilter filterId, float freqCutoff_Hz, int64_t sampleRate);

template <typename TSample>
double calcAbsMax(const TWaveformT<TSample> & waveform) {
    double amax = 0.0f;
    for (auto i = 0; i < (int) waveform.size(); ++i) if (std::abs(waveform[i]) > amax) amax = std::abs(waveform[i]);

    return amax;
}

template double calcAbsMax<TSampleF>(const TWaveformT<TSampleF> & waveform);

template <typename TSample>
bool saveToFile(const std::string & fname, TWaveformT<TSample> & waveform) {
    static_assert(std::is_same<TSample, TSampleF>::value, "Sample type not supported");

    std::ofstream fout(fname, std::ios::binary);
    if (fout.good() == false) {
        return false;
    }

    auto totalSize_bytes = sizeof(TSample)*waveform.size();

    fout.write((char *)(waveform.data()), totalSize_bytes);
    fout.close();

    printf("Total data saved: %g MB\n", ((float)(totalSize_bytes)/1024.0f/1024.0f));

    return true;
}

template bool saveToFile<TSampleF>(const std::string & fname, TWaveformT<TSampleF> & waveform);

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
template bool readFromFile<TSampleF, TSampleF>(const std::string & fname, TWaveformT<TSampleF> & res);

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
// filters
//

TFilterCoefficients calculateCoefficientsFirstOrderHighPass(int fc, int fs) {
    TFilterCoefficients res;

    float th = 2.0 * pi * fc / fs;
    float g = cos(th) / (1.0 + sin(th));
    res.a0 = (1.0 + g) / 2.0;
    res.a1 = -((1.0 + g) / 2.0);
    res.a2 = 0.0;
    res.b1 = -g;
    res.b2 = 0.0;

    return res;
}

TFilterCoefficients calculateCoefficientsSecondOrderButterworthHighPass(int fc, int fs) {
    TFilterCoefficients res;

    float c = tan(pi*fc / fs);
    res.a0 = 1.0 / (1.0 + sqrt2*c + pow(c, 2.0));
    res.a1 = -2.0 * res.a0;
    res.a2 = res.a0;
    res.b1 = 2.0 * res.a0*(pow(c, 2.0) - 1.0);
    res.b2 = res.a0 * (1.0 - sqrt2*c + pow(c, 2.0));

    return res;
}

TSampleF filterFirstOrderHighPass(TFilterCoefficients & coefficients, TSampleF sample) {
    TSampleF xn = sample;
    TSampleF yn =
        coefficients.a0*xn + coefficients.a1*coefficients.xnz1 + coefficients.a2*coefficients.xnz2 -
                             coefficients.b1*coefficients.ynz1 - coefficients.b2*coefficients.ynz2;

    coefficients.xnz2 = coefficients.xnz1;
    coefficients.xnz1 = xn;
    coefficients.ynz2 = coefficients.ynz1;
    coefficients.ynz1 = yn;

    return yn;
}

TSampleF filterSecondOrderButterworthHighPass(TFilterCoefficients & coefficients, TSampleF sample) {
    TSampleF xn = sample;
    TSampleF yn =
        coefficients.a0*xn + coefficients.a1*coefficients.xnz1 + coefficients.a2*coefficients.xnz2 -
                             coefficients.b1*coefficients.ynz1 - coefficients.b2*coefficients.ynz2;

    coefficients.xnz2 = coefficients.xnz1;
    coefficients.xnz1 = xn;
    coefficients.ynz2 = coefficients.ynz1;
    coefficients.ynz1 = yn;

    return yn;
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

template<>
std::tuple<int64_t, int64_t> calcSum(const TWaveformViewT<TSampleMI16> & waveform) {
    int64_t sum = 0;
    int64_t sum2 = 0;

    auto samples = waveform.samples;
    auto n       = waveform.n;

    for (int is = 0; is < n; ++is) {
        for (int j = 0; j < TSampleMI16::N; j++) {
            int32_t a0 = samples[is][j];
            sum += a0;
            sum2 += a0*a0;
        }
    }

    return std::tuple<int64_t, int64_t>(sum, sum2);
}

template<typename T>
std::tuple<int64_t, int64_t> calcSum(const TWaveformViewT<T> & waveform) {
    int64_t sum = 0;
    int64_t sum2 = 0;

    auto samples = waveform.samples;
    auto n       = waveform.n;

    for (int is = 0; is < n; ++is) {
        int32_t a0 = samples[is];
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

template<>
TValueCC calcCC(
    const TWaveformViewT<TSampleMI16> & waveform0,
    const TWaveformViewT<TSampleMI16> & waveform1,
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
        for (int j = 0; j < TSampleMI16::N; j++) {
            int32_t a0 = samples0[is][j];
            int32_t a1 = samples1[is][j];

            sum1 += a1;
            sum12 += a1*a1;
            sum01 += a0*a1;
        }
    }

    n *= TSampleMI16::N;

    {
        double nom   = sum01*n - sum0*sum1;
        double den2a = sum02*n - sum0*sum0;
        double den2b = sum12*n - sum1*sum1;
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
        int32_t a0 = samples0[is];
        int32_t a1 = samples1[is];

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
    int nWorkers = std::min(4, std::max(1, int(std::thread::hardware_concurrency()) - 2));
#else
    int nWorkers = std::min(4u, std::thread::hardware_concurrency());
    std::mutex mutex;
    std::vector<std::thread> workers(nWorkers);
    for (int i = 0; i < (int) workers.size(); ++i) {
        auto & worker = workers[i];
        worker = std::thread([&, sum0 = sum0, sum02 = sum02, i]() {
            TOffset cbesto = -1;
            TValueCC cbestcc = -1.0f;

            for (int o = -alignWindow + i; o <= alignWindow; o += nWorkers) {
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

    for (int o = 0; o <= 2*alignWindow; ++o) {
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

template<>
bool calculateSimilartyMap(
        const int32_t keyPressWidth_samples,
        const int32_t alignWindow_samples,
        const int32_t offsetFromPeak_samples,
        TKeyPressCollectionT<TSampleMI16> & keyPresses,
        TSimilarityMap & res) {
    int nPresses = keyPresses.size();

    int w = keyPressWidth_samples;
    int a = alignWindow_samples;

    res.clear();
    res.resize(nPresses);
    for (auto & x : res) x.resize(nPresses);

#ifdef __EMSCRIPTEN__
    int nWorkers = std::min(kMaxThreads, std::max(1, int(std::thread::hardware_concurrency()) - 2));
#else
    int nWorkers = std::thread::hardware_concurrency();
#endif

    std::vector<std::thread> workers(nWorkers);
    for (int iw = 0; iw < (int) workers.size(); ++iw) {
        auto & worker = workers[iw];
        worker = std::thread([&](int ith) {
            for (int i = ith; i < nPresses; i += nWorkers) {
                res[i][i].cc = 1.0f;
                res[i][i].offset = 0;

                const auto & waveform0 = keyPresses[i].waveform;
                const auto & pos0      = keyPresses[i].pos;

                auto & avgcc = keyPresses[i].ccAvg;

                const auto samples0 = waveform0.samples;

                for (int j = i + 1; j < nPresses; ++j) {
                    if (i == j) continue;

                    const auto waveform1 = keyPresses[j].waveform;
                    const auto pos1      = keyPresses[j].pos;

                    const auto samples1 = waveform1.samples;
                    const auto ret = findBestCC(TWaveformViewT<TSampleMI16> { samples0 + pos0 + offsetFromPeak_samples - w,     2*w },
                                                TWaveformViewT<TSampleMI16> { samples1 + pos1 + offsetFromPeak_samples - w - a, 2*w + 2*a }, a);

                    const auto bestcc     = std::get<0>(ret);
                    const auto bestoffset = std::get<1>(ret);

                    res[i][j].cc = bestcc;
                    res[i][j].offset = bestoffset;

                    res[j][i].cc = bestcc;
                    res[j][i].offset = -bestoffset;

                    avgcc += bestcc;
                }
                avgcc /= (nPresses - 1);
            }
        }, iw);
    }

    for (auto & worker : workers) worker.join();

    return true;
}

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

#ifdef __EMSCRIPTEN__
    int nWorkers = std::min(kMaxThreads, std::max(1, int(std::thread::hardware_concurrency()) - 2));
#else
    int nWorkers = std::thread::hardware_concurrency();
#endif

    std::vector<std::thread> workers(nWorkers);
    for (int iw = 0; iw < (int) workers.size(); ++iw) {
        auto & worker = workers[iw];
        worker = std::thread([&](int ith) {
            for (int i = ith; i < nPresses; i += nWorkers) {
                res[i][i].cc = 1.0f;
                res[i][i].offset = 0;

                const auto & waveform0 = keyPresses[i].waveform;
                const auto & pos0      = keyPresses[i].pos;

                auto & avgcc = keyPresses[i].ccAvg;

                const auto samples0 = waveform0.samples;

                for (int j = i + 1; j < nPresses; ++j) {
                    if (i == j) continue;

                    const auto waveform1 = keyPresses[j].waveform;
                    const auto pos1      = keyPresses[j].pos;

                    const auto samples1 = waveform1.samples;
                    const auto ret = findBestCC(TWaveformViewT<T> { samples0 + pos0 + offsetFromPeak_samples - w,     2*w },
                                                TWaveformViewT<T> { samples1 + pos1 + offsetFromPeak_samples - w - a, 2*w + 2*a }, a);

                    const auto bestcc     = std::get<0>(ret);
                    const auto bestoffset = std::get<1>(ret);

                    res[i][j].cc = bestcc;
                    res[i][j].offset = bestoffset;

                    res[j][i].cc = bestcc;
                    res[j][i].offset = -bestoffset;

                    avgcc += bestcc;
                }
                avgcc /= (nPresses - 1);
            }
        }, iw);
    }

    for (auto & worker : workers) worker.join();

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

template<>
bool findKeyPresses(
        const TWaveformViewT<TSampleMI16> & waveform,
        TKeyPressCollectionT<TSampleMI16> & res,
        TWaveformT<TSampleMI16> & waveformThreshold,
        TWaveformT<TSampleMI16> & waveformMax,
        double thresholdBackground,
        int historySize,
        int historySizeReset,
        bool removeLowPower) {
    res.clear();
    waveformThreshold.resize(waveform.n);
    waveformMax.resize(waveform.n);

    int rbBegin = 0;
    double rbAverage = 0.0;
    std::vector<double> rbSamples(8*historySize, 0.0);

    int k = historySize;
    std::deque<int64_t> que(k);

    auto samples = waveform.samples;
    auto n       = waveform.n;

    TWaveformT<TSampleMI16> waveformAbs(n);
    for (int64_t i = 0; i < n; ++i) {
        waveformAbs[i][0] = std::abs(samples[i][0]);
    }

    for (int64_t i = 0; i < n; ++i) {
        {
            int64_t ii = i - k/2;
            if (ii >= 0) {
                rbAverage *= rbSamples.size();
                rbAverage -= rbSamples[rbBegin];
                double acur = waveformAbs[i][0];
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
                double acur = waveformAbs[itest][0];
                if (acur > thresholdBackground*rbAverage) {
                    res.emplace_back(TKeyPressDataT<TSampleMI16> { std::move(waveform), itest, 0.0, -1, -1, '?' });
                }
            }
            waveformThreshold[itest][0] = thresholdBackground*rbAverage;
            waveformMax[itest] = waveformAbs[que.front()];
        }
    }

    if (removeLowPower) {
        while (true) {
            auto oldn = res.size();

            double avgPower = 0.0;
            for (const auto & kp : res) {
                avgPower += waveformAbs[kp.pos][0];
            }
            avgPower /= res.size();

            auto tmp = std::move(res);
            for (const auto & kp : tmp) {
                if (waveformAbs[kp.pos][0] > 0.3*avgPower) {
                    res.push_back(kp);
                }
            }

            if (res.size() == oldn) break;
        }
    }

    if (res.size() > 1) {
        TKeyPressCollectionT<TSampleMI16> res2;
        res2.push_back(res.front());

        for (int i = 1; i < (int) res.size(); ++i) {
            if (res[i].pos - res2.back().pos > historySizeReset || waveformMax[res[i].pos] > waveformMax[res2.back().pos]) {
                res2.push_back(res[i]);
            }
        }

        std::swap(res, res2);
    }

    return true;
}

template<typename T>
bool findKeyPresses(
        const TWaveformViewT<T> & waveform,
        TKeyPressCollectionT<T> & res,
        TWaveformT<T> & waveformThreshold,
        TWaveformT<T> & waveformMax,
        double thresholdBackground,
        int historySize,
        int historySizeReset,
        bool removeLowPower) {
    res.clear();
    waveformThreshold.resize(waveform.n);
    waveformMax.resize(waveform.n);

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
                if (acur > thresholdBackground*rbAverage) {
                    res.emplace_back(TKeyPressDataT<T> { std::move(waveform), itest, 0.0, -1, -1, '?' });
                }
            }
            waveformThreshold[itest] = thresholdBackground*rbAverage;
            waveformMax[itest] = waveformAbs[que.front()];
        }
    }

    if (removeLowPower) {
        while (true) {
            auto oldn = res.size();

            double avgPower = 0.0;
            for (const auto & kp : res) {
                avgPower += waveformAbs[kp.pos];
            }
            avgPower /= res.size();

            auto tmp = std::move(res);
            for (const auto & kp : tmp) {
                if (waveformAbs[kp.pos] > 0.3*avgPower) {
                    res.push_back(kp);
                }
            }

            if (res.size() == oldn) break;
        }
    }

    if (res.size() > 1) {
        TKeyPressCollectionT<T> res2;
        res2.push_back(res.front());

        for (int i = 1; i < (int) res.size(); ++i) {
            if (res[i].pos - res2.back().pos > historySizeReset || waveformMax[res[i].pos] > waveformMax[res2.back().pos]) {
                res2.push_back(res[i]);
            }
        }

        std::swap(res, res2);
    }

    return true;
}

template bool findKeyPresses<TSampleI16>(
        const TWaveformViewT<TSampleI16> & waveform,
        TKeyPressCollectionT<TSampleI16> & res,
        TWaveformT<TSampleI16> & waveformThreshold,
        TWaveformT<TSampleI16> & waveformMax,
        double thresholdBackground,
        int historySize,
        int historySizeReset,
        bool removeLowPower);

template<typename T>
bool saveKeyPresses(const std::string & fname, const TKeyPressCollectionT<T> & keyPresses) {
    std::ofstream fout(fname, std::ios::binary);
    int n = keyPresses.size();
    fout.write((char *)(&n), sizeof(n));
    for (int i = 0; i < n; ++i) {
        fout.write((char *)(&keyPresses[i].pos), sizeof(keyPresses[i].pos));
    }
    fout.close();

    return true;
}

template bool saveKeyPresses<TSampleI16>(const std::string & fname, const TKeyPressCollectionT<TSampleI16> & keyPresses);

template<typename T>
bool loadKeyPresses(const std::string & fname, const TWaveformViewT<T> & waveform, TKeyPressCollectionT<T> & keyPresses) {
    keyPresses.clear();

    std::ifstream fin(fname, std::ios::binary);
    int n = 0;
    fin.read((char *)(&n), sizeof(n));
    keyPresses.resize(n);
    for (int i = 0; i < n; ++i) {
        keyPresses[i].waveform = waveform;
        fin.read((char *)(&keyPresses[i].pos), sizeof(keyPresses[i].pos));
    }
    fin.close();

    return true;
}

template bool loadKeyPresses<TSampleI16>(const std::string & fname, const TWaveformViewT<TSampleI16> & waveform, TKeyPressCollectionT<TSampleI16> & keyPresses);

template<typename T>
bool dumpKeyPresses(const std::string & fname, const TKeyPressCollectionT<T> & data) {
    std::ofstream fout(fname);
    for (auto & k : data) {
        fout << k.pos << " 1" << std::endl;
    }
    fout.close();
    return true;
}

template bool dumpKeyPresses<TSampleI16>(const std::string & fname, const TKeyPressCollectionT<TSampleI16> & data);

template<typename TSample>
void cbPlayback(void * userData, uint8_t * stream, int len) {
    TPlaybackDataT<TSample> * data = (TPlaybackDataT<TSample> *)(userData);
    if (data->playing == false) {
        int offset = 0;
        TSample a = 0;
        while (len > 0) {
            memcpy(stream + offset*sizeof(a), &a, sizeof(a));
            len -= sizeof(a);
            ++offset;
        }
        return;
    }
    auto end = std::min(data->idx + TPlaybackDataT<TSample>::kSamples/data->slowDown, data->waveform.n);
    auto idx = data->idx;
    auto sidx = 0;
    for (; idx < end; ++idx) {
        TSample a = data->waveform.samples[idx];
        memcpy(stream + (sidx)*sizeof(a), &a, sizeof(a));
        len -= sizeof(a);
        ++sidx;

        if (data->slowDown == 2) {
            TSample a2 = data->waveform.samples[idx + 1];
            a = 0.5*a + 0.5*a2;
            memcpy(stream + (sidx)*sizeof(a), &a, sizeof(a));
            len -= sizeof(a);
            ++sidx;
        }
    }
    while (len > 0) {
        TSample a = 0;
        memcpy(stream + (idx - data->idx)*sizeof(a), &a, sizeof(a));
        len -= sizeof(a);
        ++idx;
    }
    data->idx = idx;
}

template void cbPlayback<TSampleI16>(void * userData, uint8_t * stream, int len);

template<typename T>
bool generateLowResWaveform(const TWaveformViewT<T> & waveform, TWaveformT<T> & waveformLowRes, int nWindow) {
    waveformLowRes.resize(waveform.n);

    int k = nWindow;
    std::deque<int64_t> que(k);

    //auto [samples, n] = waveform;
    auto samples = waveform.samples;
    auto n       = waveform.n;

    TWaveformT<T> waveformAbs(n);
    for (int64_t i = 0; i < n; ++i) {
        waveformAbs[i] = std::abs(samples[i]);
    }

    for (int64_t i = 0; i < n; ++i) {
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
            waveformLowRes[itest] = waveformAbs[que.front()];
        }
    }

    return true;
}

template bool generateLowResWaveform<TSampleI16>(const TWaveformViewT<TSampleI16> & waveform, TWaveformT<TSampleI16> & waveformLowRes, int nWindow);

template<typename T>
bool adjustKeyPresses(TKeyPressCollectionT<T> & keyPresses, TSimilarityMap & sim) {
    struct Pair {
        int i = -1;
        int j = -1;
        TValueCC cc = -1.0;

        bool operator < (const Pair & a) const { return cc > a.cc; }
    };

    bool res = false;

    int n = keyPresses.size();

    std::vector<Pair> ccpairs;
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            ccpairs.emplace_back(Pair{i, j, sim[i][j].cc});
        }
    }

    int nused = 0;
    std::vector<bool> used(n, false);

    std::sort(ccpairs.begin(), ccpairs.end());

    int npairs = ccpairs.size();
    for (int ip = 0; ip < npairs; ++ip) {
        auto & curpair = ccpairs[ip];
        int k0 = curpair.i;
        int k1 = curpair.j;
        if (used[k0] && used[k1]) continue;

        if (sim[k0][k1].offset != 0) res = true;

        if (used[k1] == false) {
            keyPresses[k1].pos += sim[k0][k1].offset;
        } else {
            keyPresses[k0].pos -= sim[k0][k1].offset;
        }

        sim[k0][k1].offset = 0;
        sim[k1][k0].offset = 0;

        if (used[k0] == false) { used[k0] = true; ++nused; }
        if (used[k1] == false) { used[k1] = true; ++nused; }

        if (nused == n) break;
    }

    return res;
}

template bool adjustKeyPresses<TSampleI16>(TKeyPressCollectionT<TSampleI16> & keyPresses, TSimilarityMap & sim);


template<typename T>
bool removeLowSimilarityKeys(TKeyPressCollectionT<T> & keyPresses, TSimilarityMap & sim, double threshold) {
    const int n = keyPresses.size();
    if (n != (int) sim.size()) {
        fprintf(stderr, "removeLowSimilarityKeys: n != sim.size()\n");
        return false;
    }

    std::vector<bool> used(n, false);

    for (int i = 0; i < n; ++i) {
        if (used[i]) continue;

        for (int j = 0 ; j < n; ++j) {
            if (i == j) continue;

            if (sim[i][j].cc > threshold) {
                used[i] = true;
                used[j] = true;
                break;
            }
        }
    }

    int nRemoved = 0;
    for (int i = 0; i < n; ++i) {
        if (used[i]) continue;

        //keyPresses.erase(keyPresses.begin() + i - nRemoved);
        ++nRemoved;
    }

    if (nRemoved == 0) {
        return true;
    }

    auto keyPresses0 = keyPresses;
    keyPresses.clear();

    for (int i = 0; i < n; ++i) {
        if (used[i]) {
            keyPresses.push_back(keyPresses0[i]);
        }
    }

    for (int i = 0; i < n; ++i) {
        auto cur = sim[i];
        sim[i].clear();
        for (int j = 0; j < n; ++j) {
            if (used[j]) {
                sim[i].push_back(cur[j]);
            }
        }
    }

    auto sim0 = sim;
    sim.clear();
    for (int i = 0; i < n; ++i) {
        if (used[i]) {
            sim.push_back(sim0[i]);
        }
    }

    return true;
}

template bool removeLowSimilarityKeys<TSampleI16>(TKeyPressCollectionT<TSampleI16> & keyPresses, TSimilarityMap & sim, double threshold);
