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
#include <algorithm>

#define MY_DEBUG

using TSum                  = int64_t;
using TSum2                 = int64_t;
using TCC                   = double;
using TOffset               = int64_t;
using TMatch                = std::tuple<TCC, TOffset>;
using TSimilarityMap        = std::vector<std::vector<TMatch>>;

using TClusterId            = int32_t;
using TSampleInput          = float;
using TSample               = int32_t;
using TWaveform             = std::vector<TSample>;
using TWaveformView         = std::tuple<const TSample *, int64_t>;
using TKeyPressPosition     = int64_t;
using TKeyPressData         = std::tuple<TWaveformView, TKeyPressPosition, TClusterId, TCC>;
using TKeyPressCollection   = std::vector<TKeyPressData>;

template <typename T>
float toSeconds(T t0, T t1) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()/1024.0f;
}

inline float frand() { return ((float)rand())/RAND_MAX; }

TWaveformView getView(const TWaveform & waveform, int64_t idx) { return TWaveformView { waveform.data() + idx, waveform.size() - idx }; }

TWaveformView getView(const TWaveform & waveform, int64_t idx, int64_t len) { return TWaveformView { waveform.data() + idx, len }; }

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
    std::array<double, 4*1024> rbSamples;
    rbSamples.fill(0.0);

    int k = 1024;
    double thresholdBackground = 10.0;

    std::deque<int64_t> que(k);
    //auto [samples, nSamples] = waveform;
    auto samples  = std::get<0>(waveform);
    auto nSamples = std::get<1>(waveform);

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
                    res.push_back(TKeyPressData {waveform, itest, 0, 0.0});
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
        //auto [_k0, idx, cid, _k3] = k;
        auto idx = std::get<1>(k);
        fout << idx << " 1" << std::endl;
    }
    fout.close();
    return true;
}

std::tuple<TSum, TSum2> calcSum(const TWaveformView & waveform) {
    TSum sum = 0.0f;
    TSum2 sum2 = 0.0f;
    //auto [samples, n] = waveform;
    auto samples = std::get<0>(waveform);
    auto n       = std::get<1>(waveform);

    for (int64_t is = 0; is < n; ++is) {
        auto a0 = samples[is];
        sum += a0;
        sum2 += a0*a0;
    }

    return std::tuple<TSum, TSum2>(sum, sum2);
}

TCC calcCC(
    const TWaveformView & waveform0,
    const TWaveformView & waveform1,
    TSum sum0, TSum2 sum02) {
    TCC cc = -1.0f;

    TSum sum1 = 0.0f;
    TSum2 sum12 = 0.0f;
    TSum2 sum01 = 0.0f;

    //auto [samples0, n0] = waveform0;
    auto samples0 = std::get<0>(waveform0);
    auto n0       = std::get<1>(waveform0);

    //auto [samples1, n1] = waveform1;
    auto samples1 = std::get<0>(waveform1);
    auto n1       = std::get<1>(waveform1);

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

    //auto [samples0, n0] = waveform0;
    //auto samples0 = std::get<0>(waveform0);
    auto n0       = std::get<1>(waveform0);

    //auto [samples1, n1] = waveform1;
    auto samples1 = std::get<0>(waveform1);
    auto n1       = std::get<1>(waveform1);

#ifdef MY_DEBUG
    if (n0 + 2*alignWindow != n1) {
        printf("BUG 924830jm92, n0 = %d, n1 = %d, a = %d\n", (int) n0, (int) n1, (int) alignWindow);
    }
#endif

    //auto [sum0, sum02] = calcSum(waveform0);
    auto ret = calcSum(waveform0);
    auto sum0  = std::get<0>(ret);
    auto sum02 = std::get<1>(ret);

    for (int o = 0; o < 2*alignWindow; ++o) {
        auto cc = calcCC(waveform0, TWaveformView { samples1 + o, n0 }, sum0, sum02);
        if (cc > bestcc) {
            besto = o - alignWindow;
            bestcc = cc;
        }
    }

    return std::tuple<TCC, TOffset>(bestcc, besto);
}

bool calculateSimilartyMap(TKeyPressCollection & keyPresses, TSimilarityMap & res) {
    res.clear();
    int nPresses = keyPresses.size();

    int w = 256;
    int alignWindow = 256;

    res.resize(nPresses);
    for (auto & x : res) x.resize(nPresses);

    for (int i = 0; i < nPresses; ++i) {
        res[i][i] = TMatch { 1.0f, 0 };
        //auto & [waveform0, pos0, _i2, avgcc] = keyPresses[i];
        auto & waveform0 = std::get<0>(keyPresses[i]);
        auto & pos0      = std::get<1>(keyPresses[i]);
        //   &           = std::get<2>(keyPresses[i]);
        auto & avgcc     = std::get<3>(keyPresses[i]);

        //auto [samples0, n0] = waveform0;
        auto samples0 = std::get<0>(waveform0);
        //auto n0       = std::get<1>(waveform0);

        for (int j = 0; j < nPresses; ++j) {
            if (i == j) continue;

            //auto [waveform1, pos1, _j2, _j3] = keyPresses[j];
            auto waveform1 = std::get<0>(keyPresses[j]);
            auto pos1      = std::get<1>(keyPresses[j]);

            //auto [samples1, n1] = waveform1;
            auto samples1 = std::get<0>(waveform1);
            //auto n1       = std::get<1>(waveform1);

            //auto [bestcc, bestoffset] = findBestCC({ samples0 + pos0 + (int)(0.5f*w),               2*w },
            //                                       { samples1 + pos1 + (int)(0.5f*w) - alignWindow, 2*w + 2*alignWindow }, alignWindow);
            auto ret = findBestCC(TWaveformView { samples0 + pos0 + (int)(0.5f*w),               2*w },
                                  TWaveformView { samples1 + pos1 + (int)(0.5f*w) - alignWindow, 2*w + 2*alignWindow }, alignWindow);
            auto bestcc     = std::get<0>(ret);
            auto bestoffset = std::get<1>(ret);

            res[j][i] = TMatch { bestcc, bestoffset };
            //res[i][j] = { bestcc, -bestoffset };

            avgcc += bestcc;
        }
        avgcc /= (nPresses - 1);
    }

    return true;
}

bool clusterDBSCAN(const TSimilarityMap & sim, TCC epsCC, int minPts, TKeyPressCollection & keyPresses) {
    int n = keyPresses.size();
    for (int i = 0; i < n; ++i) {
        //auto & [_i0, _i1, cid, _i3] = keyPresses[i];
        auto & cid = std::get<2>(keyPresses[i]);

        cid = -1;
    }

    int curId = 0;

    for (int i = 0; i < n; ++i) {
        //auto & [_i0, i1, cid, _i2] = keyPresses[i];
        auto & cid = std::get<2>(keyPresses[i]);

        if (cid != -1) continue;
        std::vector<int> nbi;
        for (int j = 0; j < n; ++j) {
            //auto & [cc, x0] = sim[i][j];
            auto & cc = std::get<0>(sim[i][j]);

            if (cc > epsCC) nbi.push_back(j);
        }

        if (nbi.size() < minPts) {
            cid = 0;
            continue;
        }

        cid = ++curId;
        for (int q = 0; q < (int) nbi.size(); ++q) {
            //auto & [_q0, _q1, qcid, _q3] = keyPresses[nbi[q]];
            auto & qcid = std::get<2>(keyPresses[nbi[q]]);

            if (qcid == 0) qcid = curId;
            if (qcid != -1) continue;
            qcid = curId;

            std::vector<int> nbq;
            for (int j = 0; j < n; ++j) {
                //auto & [cc, x0] = sim[nbi[q]][j];
                auto & cc = std::get<0>(sim[nbi[q]][j]);

                if (cc > epsCC) nbq.push_back(j);
            }

            if (nbq.size() >= minPts) {
                nbi.insert(nbi.end(), nbq.begin(), nbq.end());
            }
        }
    }

    return true;
}

bool clusterSimple(const TSimilarityMap & sim, TCC tholdCC, TKeyPressCollection & keyPresses) {
    int n = keyPresses.size();
    for (int i = 0; i < n; ++i) {
        //auto & [_i0, _i1, cid, _i2] = keyPresses[i];
        auto & cid = std::get<2>(keyPresses[i]);

        cid = 0;
    }

    int curId = 0;

    int nIter = 1000000;
    for (int iter = 0; iter < nIter; ++iter) {
        int i = rand()%n;
        int j = rand()%n;
        while (i == j) {
            i = rand()%n;
            j = rand()%n;
        }

        //auto & [_i0, i1, icid, _i3] = keyPresses[i];
        auto & icid = std::get<2>(keyPresses[i]);

        //auto & [_j0, j1, jcid, _j3] = keyPresses[j];
        auto & jcid = std::get<2>(keyPresses[j]);

        //if (icid != 0 && jcid != 0) continue;

        //auto & [cc, x] = sim[i][j];
        auto & cc = std::get<0>(sim[i][j]);

        if (cc < tholdCC) continue;
        auto r = frand();
        if (std::pow(r, 1.0) < cc) {
            if (icid == jcid && icid == 0) {
                icid = jcid = ++curId;
            } else if (icid == 0) {
                icid = jcid;
            } else {
                jcid = icid;
            }
        }

        //if (icid == jcid) {
        //    //if (std::pow(r, 1.0) > cc) {
        //    //    icid = i + 1;
        //    //    jcid = j + 1;
        //    //}
        //} else {
        //    if (std::pow(r, 1.0) < cc) {
        //        icid = jcid;
        //    }
        //}
    }

    for (int i = 0; i < n; ++i) {
        //auto & [_i0, _i1, icid, _i3] = keyPresses[i];
        auto & icid = std::get<2>(keyPresses[i]);

        if (icid == 0) icid = ++curId;
        icid = -icid;
    }

    curId = 0;
    for (int i = 0; i < n; ++i) {
        //auto [_i0, i1, icid, _i3] = keyPresses[i];
        auto & icid = std::get<2>(keyPresses[i]);

        if (icid > 0) continue;
        for (int j = i; j < n; ++j) {
            //auto & [_j0, _j1, jcid, _j3] = keyPresses[j];
            auto & jcid = std::get<2>(keyPresses[j]);

            if (jcid == icid) {
                jcid = curId;
            }
        }
        ++curId;
    }

    return true;
}

bool clusterG(const TSimilarityMap & sim, TKeyPressCollection & keyPresses) {
    struct Pair {
        int i = -1;
        int j = -1;
        TCC cc = -1.0;

        bool operator < (const Pair & a) const { return cc > a.cc; }
    };

    int n = keyPresses.size();

    int nclusters = 0;
    for (int i = 0; i < n; ++i) {
        //auto & [_i0, _i1, ci, _i3] = keyPresses[i];
        auto & ci = std::get<2>(keyPresses[i]);

        ci = i + 1;
        ++nclusters;
    }

    std::vector<Pair> ccpairs;
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            //auto & [cc, _] = sim[i][j];
            auto & cc = std::get<0>(sim[i][j]);

            ccpairs.emplace_back(Pair{i, j, cc});
        }
    }

    std::sort(ccpairs.begin(), ccpairs.end());

    printf("[+] Top 10 pairs\n");
    for (int i = 0; i < 10; ++i) {
        printf("    Pair %d: %d %d %g\n", i, ccpairs[i].i, ccpairs[i].j, ccpairs[i].cc);
    }

    int npairs = ccpairs.size();
    for (int ip = 0; ip < npairs; ++ip) {
        auto & curpair = ccpairs[ip];
        //if (frand() > curpair.cc) continue;
        if (curpair.cc < 0.50) break;

        //auto [_i0, __i0, ci, _i3] = keyPresses[curpair.i];
        auto ci = std::get<2>(keyPresses[curpair.i]);

        //auto [_j1, __j1, cj, _j3] = keyPresses[curpair.j];
        auto cj = std::get<2>(keyPresses[curpair.j]);

        if (ci == cj) continue;
        auto cnew = std::min(ci, cj);
        int nsum = 0;
        int nsumi = 0;
        int nsumj = 0;
        double sumcc = 0.0;
        double sumcci = 0.0;
        double sumccj = 0.0;
        for (int k = 0; k < n; ++k) {
            //auto & [_k0, _k1, ck, _k3] = keyPresses[k];
            auto ck = std::get<2>(keyPresses[k]);

            for (int q = 0; q < n; ++q) {
                if (q == k) continue;
                //auto & [_q0, _q1, cq, _q3] = keyPresses[q];
                auto cq = std::get<2>(keyPresses[q]);

                if ((ck == ci || ck == cj) && (cq == ci || cq == cj)) {
                    //auto & [cc, offset] = sim[k][q];
                    auto & cc     = std::get<0>(sim[k][q]);
                    //auto & offset = std::get<1>(sim[k][q]);

                    sumcc += cc;
                    ++nsum;
                }
                if (ck == ci && cq == ci) {
                    //auto & [cc, offset] = sim[k][q];
                    auto & cc     = std::get<0>(sim[k][q]);
                    //auto & offset = std::get<1>(sim[k][q]);

                    sumcci += cc;
                    ++nsumi;
                }
                if (ck == cj && cq == cj) {
                    //auto & [cc, offset] = sim[k][q];
                    auto & cc     = std::get<0>(sim[k][q]);
                    //auto & offset = std::get<1>(sim[k][q]);

                    sumccj += cc;
                    ++nsumj;
                }
            }
        }
        sumcc /= nsum;
        if (nsumi > 0) sumcci /= nsumi;
        if (nsumj > 0) sumccj /= nsumj;
        printf("Merge avg n = %4d, cc = %8.5f, ni = %4d, cci = %8.5f, nj = %4d, ccj = %8.5f\n", nsum, sumcc, nsumi, sumcci, nsumj, sumccj);

        //if (sumcc < 1.000*curpair.cc) continue;
        //if (sumcc > 0.75*sumccj && sumcc > 0.75*sumcci) {
        if (sumcc > 0.4*(sumcci + sumccj)) {
        } else {
            continue;
        }

        for (int k = 0; k < n; ++k) {
            //auto & [_k0, _k1, ck, _k3] = keyPresses[k];
            auto & ck = std::get<2>(keyPresses[k]);

            if (ck == ci || ck == cj) ck = cnew;
        }
        --nclusters;

        printf("Clusters %3d %5.4f:", nclusters, curpair.cc);
        for (int k = 0; k < n; ++k) {
            //auto & [_k0, _k1, ck, _k3] = keyPresses[k];
            auto & ck = std::get<2>(keyPresses[k]);

            printf(" %3d", ck);
        }
        printf("\n");
    }

    return true;
}

bool clusterG2(TSimilarityMap & sim, TKeyPressCollection & keyPresses) {
    struct Pair {
        int i = -1;
        int j = -1;
        TCC cc = -1.0;

        bool operator < (const Pair & a) const { return cc > a.cc; }
    };

    int n = keyPresses.size();

    int nclusters = 0;
    for (int i = 0; i < n; ++i) {
        //auto & [_i0, _i1, ci, _i3] = keyPresses[i];
        auto & ci = std::get<2>(keyPresses[i]);

        ci = i + 1;
        ++nclusters;
    }

    std::vector<double> nsum(n + 1, 1.0);
    std::vector<Pair> ccpairs;
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            //auto & [cc, _] = sim[i][j];
            auto & cc = std::get<0>(sim[i][j]);

            ccpairs.emplace_back(Pair{i, j, cc});
        }
    }

    std::sort(ccpairs.begin(), ccpairs.end());

    printf("[+] Top 10 pairs\n");
    for (int i = 0; i < 10; ++i) {
        printf("    Pair %d: %d %d %g\n", i, ccpairs[i].i, ccpairs[i].j, ccpairs[i].cc);
    }

    int npairs = ccpairs.size();
    for (int ip = 0; ip < npairs; ++ip) {
        auto & curpair = ccpairs[ip];
        if (frand() > curpair.cc) continue;
        if (curpair.cc < 0.50) break;

        //auto [_i0, __i0, ci, _i3] = keyPresses[curpair.i];
        auto & ci = std::get<2>(keyPresses[curpair.i]);

        //auto [_j1, __j1, cj, _j3] = keyPresses[curpair.j];
        auto & cj = std::get<2>(keyPresses[curpair.j]);

        if (ci == cj) continue;
        auto cnew = std::min(ci, cj);
        for (int k = 0; k < n; ++k) {
            //auto & [_k0, _k1, ck, _k3] = keyPresses[k];
            auto & ck = std::get<2>(keyPresses[k]);

            if (ck == ci || ck == cj) ck = cnew;
        }
        --nclusters;

        printf("Clusters %3d %5.4f:", nclusters, curpair.cc);
        for (int k = 0; k < n; ++k) {
            //auto & [_k0, _k1, ck, _k3] = keyPresses[k];
            auto & ck = std::get<2>(keyPresses[k]);

            printf(" %3d", ck);
        }
        printf("\n");

        if (nclusters == 56) break;
    }

    return true;
}

int main(int argc, char ** argv) {
    srand(time(0));

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

    //{
    //    std::ofstream fout("waveform.plot");
    //    for (auto & a : waveformInput) {
    //        fout << a << std::endl;
    //    }
    //    fout.close();
    //}

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
        for (auto & k : keyPresses) {
            //auto & [_k0, pos, _k2, _k3] = k;
            auto & pos = std::get<1>(k);

            printf("    position - %d\n", (int) pos);
        }
        printf("[+] Search took %4.3f seconds\n", toSeconds(tStart, tEnd));

        dumpKeyPresses("key_presses.plot", keyPresses);
    }

    //keyPresses.erase(keyPresses.begin(), keyPresses.begin() + keyPresses.size()/2);

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

    int n = keyPresses.size();
    for (int i = 0; i < n; ++i) {
        //auto [_i0, _i1, _i2, avgcc] = keyPresses[i];
        auto & avgcc = std::get<3>(keyPresses[i]);

        printf("    Average CC for keypress %4d - %6.3f\n", i, avgcc);
    }

    printf("%5d ", -1);
    for (int j = 0; j < n; ++j) {
        printf("%4d ", j);
    }
    printf("\n");
    printf("--------------------------------------------------------------------------------------------------------------------------------------\n");

    for (int i = 0; i < n; ++i) {
        printf("%2d  | ", i);
        for (int j = 0; j < n; ++j) {
            //auto [cc, offset] = similarityMap[i][j];
            auto cc     = std::get<0>(similarityMap[i][j]);
            //auto offset = std::get<1>(similarityMap[i][j]);

            if (cc > -0.45) {
                printf("%4.0f ", cc*100);
            } else {
                printf("     ");
            }
            //printf("%3d - %3d -> %4.3f, %3d\n", i, j, cc, (int) offset);
        }
        printf("\n");
    }

    printf("\n");

    //clusterDBSCAN(similarityMap, 0.8, 1, keyPresses);
    //clusterSimple(similarityMap, 0.50, keyPresses);
    clusterG(similarityMap, keyPresses);

    for (const auto & k : keyPresses) {
        //auto & [_k0, _k1, cid, _k2] = k;
        auto & cid = std::get<2>(k);

        //printf("%c", 'a' + cid - 1);
        printf("%d ", cid);
    }
    printf("\n");

    for (const auto & k : keyPresses) {
        //auto & [_k0, _k1, cid, _k3] = k;
        auto & cid = std::get<2>(k);

        printf("%c", 'a' + cid - 1);
    }
    printf("\n");

    return 0;
}
