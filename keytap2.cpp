/*! \file keytap2.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "common.h"

#include <chrono>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#define MY_DEBUG

using TSampleInput          = TSampleF;
using TSample               = TSampleI16;
using TWaveform             = TWaveformI16;
using TWaveformView         = TWaveformViewI16;
using TKeyPressData         = TKeyPressDataI16;
using TKeyPressCollection   = TKeyPressCollectionI16;

bool clusterDBSCAN(const TSimilarityMap & sim, TValueCC epsCC, int minPts, TKeyPressCollection & keyPresses) {
    int n = keyPresses.size();
    for (int i = 0; i < n; ++i) {
        auto & cid = keyPresses[i].cid;

        cid = -1;
    }

    int curId = 0;

    for (int i = 0; i < n; ++i) {
        auto & cid = keyPresses[i].cid;

        if (cid != -1) continue;
        std::vector<int> nbi;
        for (int j = 0; j < n; ++j) {
            auto & cc = sim[i][j].cc;

            if (cc > epsCC) nbi.push_back(j);
        }

        if ((int) nbi.size() < minPts) {
            cid = 0;
            continue;
        }

        cid = ++curId;
        for (int q = 0; q < (int) nbi.size(); ++q) {
            auto & qcid = keyPresses[nbi[q]].cid;

            if (qcid == 0) qcid = curId;
            if (qcid != -1) continue;
            qcid = curId;

            std::vector<int> nbq;
            for (int j = 0; j < n; ++j) {
                auto & cc = sim[nbi[q]][j].cc;

                if (cc > epsCC) nbq.push_back(j);
            }

            if ((int) nbq.size() >= minPts) {
                nbi.insert(nbi.end(), nbq.begin(), nbq.end());
            }
        }
    }

    return true;
}

bool clusterSimple(const TSimilarityMap & sim, TValueCC tholdCC, TKeyPressCollection & keyPresses) {
    int n = keyPresses.size();
    for (int i = 0; i < n; ++i) {
        auto & cid = keyPresses[i].cid;

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

        auto & icid = keyPresses[i].cid;
        auto & jcid = keyPresses[j].cid;

        //if (icid != 0 && jcid != 0) continue;

        auto & cc = sim[i][j].cc;

        if (cc < tholdCC) continue;
        auto r = frand();
        if (r < cc) {
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
        auto & icid = keyPresses[i].cid;

        if (icid == 0) icid = ++curId;
        icid = -icid;
    }

    curId = 0;
    for (int i = 0; i < n; ++i) {
        auto & icid = keyPresses[i].cid;

        if (icid > 0) continue;
        for (int j = i; j < n; ++j) {
            auto & jcid = keyPresses[j].cid;

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
        TValueCC cc = -1.0;

        bool operator < (const Pair & a) const { return cc > a.cc; }
    };

    int n = keyPresses.size();

    int nclusters = 0;
    for (int i = 0; i < n; ++i) {
        auto & ci = keyPresses[i].cid;

        ci = i + 1;
        ++nclusters;
    }

    std::vector<Pair> ccpairs;
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            auto & cc = sim[i][j].cc;

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

        auto ci = keyPresses[curpair.i].cid;
        auto cj = keyPresses[curpair.j].cid;

        if (ci == cj) continue;
        auto cnew = std::min(ci, cj);
        int nsum = 0;
        int nsumi = 0;
        int nsumj = 0;
        double sumcc = 0.0;
        double sumcci = 0.0;
        double sumccj = 0.0;
        for (int k = 0; k < n; ++k) {
            auto ck = keyPresses[k].cid;

            for (int q = 0; q < n; ++q) {
                if (q == k) continue;
                auto cq = keyPresses[q].cid;

                if ((ck == ci || ck == cj) && (cq == ci || cq == cj)) {
                    auto & cc     = sim[k][q].cc;
                    //auto & offset = sim[k][q].offset;

                    sumcc += cc;
                    ++nsum;
                }
                if (ck == ci && cq == ci) {
                    auto & cc     = sim[k][q].cc;
                    //auto & offset = sim[k][q].offset;

                    sumcci += cc;
                    ++nsumi;
                }
                if (ck == cj && cq == cj) {
                    auto & cc     = sim[k][q].cc;
                    //auto & offset = sim[k][q].offset;

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
            auto & ck = keyPresses[k].cid;

            if (ck == ci || ck == cj) ck = cnew;
        }
        --nclusters;

        printf("Clusters %3d %5.4f:", nclusters, curpair.cc);
        for (int k = 0; k < n; ++k) {
            auto & ck = keyPresses[k].cid;

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
        TValueCC cc = -1.0;

        bool operator < (const Pair & a) const { return cc > a.cc; }
    };

    int n = keyPresses.size();

    int nclusters = 0;
    for (int i = 0; i < n; ++i) {
        auto & ci = keyPresses[i].cid;

        ci = i + 1;
        ++nclusters;
    }

    std::vector<double> nsum(n + 1, 1.0);
    std::vector<Pair> ccpairs;
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            auto & cc = sim[i][j].cc;

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

        auto & ci = keyPresses[curpair.i].cid;
        auto & cj = keyPresses[curpair.j].cid;

        if (ci == cj) continue;
        auto cnew = std::min(ci, cj);
        for (int k = 0; k < n; ++k) {
            auto & ck = keyPresses[k].cid;

            if (ck == ci || ck == cj) ck = cnew;
        }
        --nclusters;

        printf("Clusters %3d %5.4f:", nclusters, curpair.cc);
        for (int k = 0; k < n; ++k) {
            auto & ck = keyPresses[k].cid;

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
    if (readFromFile<TSampleF>(argv[1], waveformInput) == false) {
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
    TWaveform waveformThreshold;
    TWaveform waveformMax;
    {
        auto tStart = std::chrono::high_resolution_clock::now();
        printf("[+] Searching for key presses\n");
        if (findKeyPresses(getView(waveformInput, 0), keyPresses, waveformThreshold, waveformMax, 10.0, 6*1024) == false) {
            printf("Failed to detect keypresses\n");
            return -2;
        }
        auto tEnd = std::chrono::high_resolution_clock::now();
        printf("[+] Detected a total of %d potential key presses\n", (int) keyPresses.size());
        for (auto & k : keyPresses) {
            auto & pos = k.pos;

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
        if (calculateSimilartyMap(256, 256, 0, keyPresses, similarityMap) == false) {
            printf("Failed to calculate similariy map\n");
            return -3;
        }
        auto tEnd = std::chrono::high_resolution_clock::now();
        printf("[+] Calculation took %4.3f seconds\n", toSeconds(tStart, tEnd));
    }

    int n = keyPresses.size();
    for (int i = 0; i < n; ++i) {
        auto & avgcc = keyPresses[i].ccAvg;

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
            auto cc     = similarityMap[i][j].cc;
            //auto offset = similarityMap[i][j].offset;

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
        auto & cid = k.cid;

        //printf("%c", 'a' + cid - 1);
        printf("%d ", cid);
    }
    printf("\n");

    for (const auto & k : keyPresses) {
        auto & cid = k.cid;

        printf("%c", 'a' + cid - 1);
    }
    printf("\n");

    return 0;
}
