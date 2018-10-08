/*! \file keytap2.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <deque>
#include <fstream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#define MY_DEBUG

using TSampleInput          = float;
using TSample               = int32_t;
using TWaveform             = std::vector<TSample>;
using TWaveformView         = std::tuple<const TSample *, int64_t>;
using TKeyPressPosition     = int64_t;
using TKeyPressData         = std::tuple<TWaveformView, TKeyPressPosition>;
using TKeyPressCollection   = std::vector<TKeyPressData>;

using TSum                  = int64_t;
using TSum2                 = int64_t;
using TCC                   = double;
using TOffset               = int64_t;
using TMatch                = std::tuple<TCC, TOffset>;
using TSimilarityMap        = std::vector<std::vector<TMatch>>;

template <typename T>
float toSeconds(T t0, T t1) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()/1024.0f;
}

TWaveformView getView(const TWaveform & waveform, int64_t idx) { return { waveform.data() + idx, waveform.size() - idx }; }

TWaveformView getView(const TWaveform & waveform, int64_t idx, int64_t len) { return { waveform.data() + idx, len }; }

bool readFromFile(const std::string & fname, TWaveform & res) {
    std::ifstream fin(fname, std::ios::binary | std::ios::ate);
    if (fin.good() == false) {
        return false;
    }

    {
        std::streamsize size = fin.tellg();
        fin.seekg(0, std::ios::beg);

        static_assert(std::is_same<TSampleInput, float>::value, "TSampleInput not recognised");
        static_assert(
            std::is_same<TSample, float>::value
            || std::is_same<TSample, int16_t>::value
            || std::is_same<TSample, int32_t>::value
                      , "TSampleInput not recognised");

        if (std::is_same<TSample, int16_t>::value) {
            std::vector<TSampleInput> buf(size/sizeof(TSampleInput));
            res.resize(size/sizeof(TSampleInput));
            fin.read((char *)(buf.data()), size);
            double amax = 0.0f;
            for (auto i = 0; i < buf.size(); ++i) if (std::abs(buf[i]) > amax) amax = std::abs(buf[i]);
            for (auto i = 0; i < buf.size(); ++i) res[i] = std::round(32000.0*(buf[i]/amax));
            //double asum = 0.0f;
            //for (auto i = 0; i < buf.size(); ++i) asum += std::abs(buf[i]); asum /= buf.size(); asum *= 10.0;
            //for (auto i = 0; i < buf.size(); ++i) res[i] = std::round(2000.0*(buf[i]/asum));
        } else if (std::is_same<TSample, int32_t>::value) {
            std::vector<TSampleInput> buf(size/sizeof(TSampleInput));
            res.resize(size/sizeof(TSampleInput));
            fin.read((char *)(buf.data()), size);
            double amax = 0.0f;
            for (auto i = 0; i < buf.size(); ++i) if (std::abs(buf[i]) > amax) amax = std::abs(buf[i]);
            for (auto i = 0; i < buf.size(); ++i) res[i] = std::round(32000.0*(buf[i]/amax));
        } else if (std::is_same<TSample, float>::value) {
            res.resize(size/sizeof(TSample));
            fin.read((char *)(res.data()), size);
        } else {
        }
    }

    fin.close();

    return true;
}

bool findKeyPresses(const TWaveformView & waveform, TKeyPressCollection & res) {
    res.clear();

    int rbBegin = 0;
    double rbAverage = 0.0;
    std::array<double, 1024> rbSamples;
    rbSamples.fill(0.0);

    int k = 512;
    double thresholdBackground = 10.0;

    std::deque<int64_t> que(k);
    auto [samples, nSamples] = waveform;
    for (int64_t i = 0; i < nSamples; ++i) {
        {
            rbAverage *= rbSamples.size();
            rbAverage -= rbSamples[rbBegin];
            double acur = std::abs(samples[i]);
            rbSamples[rbBegin] = acur;
            rbAverage += acur;
            rbAverage /= rbSamples.size();
            if (++rbBegin >= rbSamples.size()) {
                rbBegin = 0;
            }
        }

        if (i < k) {
            while((!que.empty()) && samples[i] >= samples[que.back()]) {
                que.pop_back();
            }
            que.push_back(i);
        } else {
            while((!que.empty()) && que.front() <= i - k) {
                que.pop_front();
            }

            while((!que.empty()) && samples[i] >= samples[que.back()]) {
                que.pop_back();
            }

            que.push_back(i);

            int64_t itest = i - k/2;
            if (itest >= 2*k && itest < nSamples - 2*k && que.front() == itest) {
                double acur = samples[itest];
                if (acur > thresholdBackground*rbAverage){
                    res.push_back({waveform, itest});
                }
            }
        }
    }

    return true;
}

bool findKeyPresses(const TWaveform & waveform, TKeyPressCollection & res) {
    return findKeyPresses(getView(waveform, 0), res);
}

bool dumpKeyPresses(const std::string & fname, const TKeyPressCollection & data) {
    std::ofstream fout(fname);
    for (auto & k : data) {
        auto [waveform, idx] = k;
        fout << idx << " 1" << std::endl;
    }
    fout.close();
    return true;
}

std::tuple<TSum, TSum2> calcSum(const TWaveformView & waveform) {
    TSum sum = 0.0f;
    TSum2 sum2 = 0.0f;
    auto [samples, n] = waveform;
    for (int64_t is = 0; is < n; ++is) {
        auto a0 = samples[is];
        sum += a0;
        sum2 += a0*a0;
    }

    return { sum, sum2 };
}

TCC calcCC(
    const TWaveformView & waveform0,
    const TWaveformView & waveform1,
    TSum sum0, TSum2 sum02) {
    TCC cc = -1.0f;

    TSum sum1 = 0.0f;
    TSum2 sum12 = 0.0f;
    TSum2 sum01 = 0.0f;

    auto [samples0, n0] = waveform0;
    auto [samples1, n1] = waveform1;

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

std::tuple<TCC, TOffset> findBestCC(
    const TWaveformView & waveform0,
    const TWaveformView & waveform1,
    int64_t alignWindow) {
    TCC bestcc = -1.0;
    TOffset besto = -1;

    auto [samples0, n0] = waveform0;
    auto [samples1, n1] = waveform1;

#ifdef MY_DEBUG
    if (n0 + 2*alignWindow != n1) {
        printf("BUG 924830jm92, n0 = %d, n1 = %d, a = %d\n", (int) n0, (int) n1, (int) alignWindow);
    }
#endif

    auto [sum0, sum02] = calcSum(waveform0);

    for (int o = 0; o < 2*alignWindow; ++o) {
        auto cc = calcCC(waveform0, { samples1 + o, n0 }, sum0, sum02);
        if (cc > bestcc) {
            besto = o - alignWindow;
            bestcc = cc;
        }
    }

    return { bestcc, besto };
}
bool calculateSimilartyMap(const TKeyPressCollection & keyPresses, TSimilarityMap & res) {
    res.clear();
    int nPresses = keyPresses.size();

    int w = 512;
    int alignWindow = 64;

    res.resize(nPresses);
    for (auto & x : res) x.resize(nPresses);

    for (int i = 0; i < nPresses; ++i) {
        res[i][i] = { 1.0f, 0 };
        auto [waveform0, pos0] = keyPresses[i];
        auto [samples0, n0] = waveform0;

        for (int j = i + 1; j < nPresses; ++j) {
            auto [waveform1, pos1] = keyPresses[j];
            auto [samples1, n1] = waveform1;
            auto [bestcc, bestoffset] = findBestCC({ samples0 + pos0 - w,               2*w },
                                                   { samples1 + pos1 - w - alignWindow, 2*w + 2*alignWindow }, alignWindow);

            res[j][i] = { bestcc, bestoffset };
            res[i][j] = { bestcc, -bestoffset };
        }
    }

    return true;
}

int main(int argc, char ** argv) {
    printf("Usage: %s record.kbd\n", argv[0]);
    if (argc < 2) {
        return -1;
    }

    int64_t sampleRate = 24000;

    TWaveform waveformInput;
    printf("[+] Loading recording from '%s'\n", argv[1]);
    if (readFromFile(argv[1], waveformInput) == false) {
        printf("Specified file '%s' does not exist\n", argv[1]);
        return -1;
    }

    printf("[+] Loaded recording: of %d samples (sample size = %d bytes)\n", (int) waveformInput.size(), (int) sizeof(TSample));
    printf("    Size in memory:          %g MB\n", (float)(sizeof(TSample)*waveformInput.size())/1024/1024);
    printf("    Sample size:             %d\n", (int) sizeof(TSample));
    printf("    Total number of samples: %d\n", (int) waveformInput.size());
    printf("    Recording length:        %g seconds\n", (float)(waveformInput.size())/sampleRate);

    TKeyPressCollection keyPresses;
    {
        auto tStart = std::chrono::high_resolution_clock::now();
        printf("[+] Searching for key presses\n");
        if (findKeyPresses(waveformInput, keyPresses) == false) {
            printf("Failed to detect keypresses\n");
            return -2;
        }
        auto tEnd = std::chrono::high_resolution_clock::now();
        printf("[+] Detected a total of %d potential key presses\n", (int) keyPresses.size());
        printf("[+] Search took %4.3f seconds\n", toSeconds(tStart, tEnd));

        dumpKeyPresses("key_presses.plot", keyPresses);
    }

    //keyPresses.resize(20);

    TSimilarityMap similarityMap;
    {
        auto tStart = std::chrono::high_resolution_clock::now();
        printf("[+] Calculating CC similarity map\n");
        if (calculateSimilartyMap(keyPresses, similarityMap) == false) {
            printf("Failed to calculate similariy map\n");
            return -3;
        }
        auto tEnd = std::chrono::high_resolution_clock::now();
        printf("[+] Calculation took %4.3f seconds\n", toSeconds(tStart, tEnd));
    }

    //int n = keyPresses.size();
    //for (int i = 0; i < n; ++i) {
    //    for (int j = 0; j < n; ++j) {
    //        auto [cc, offset] = similarityMap[i][j];
    //        printf("%3d - %3d -> %4.3f, %3d\n", i, j, cc, (int) offset);
    //    }
    //}

    //return 0;
}
