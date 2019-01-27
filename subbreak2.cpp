/*! \file subbreak2.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "subbreak2.h"

#include <array>
#include <cstdlib>
#include <fstream>
#include <chrono>

namespace Cipher {
    inline float frand() { return ((float)rand())/RAND_MAX; }

    static const std::array<float, 26> kEnglishLetterFreq = {
        8.167,
        1.492,
        2.782,
        4.253,
        12.702,
        2.228,
        2.015,
        6.094,
        6.966,
        0.153,
        0.772,
        4.025,
        2.406,
        6.749,
        7.507,
        1.929,
        0.095,
        5.987,
        6.327,
        9.056,
        2.758,
        0.978,
        2.360,
        0.150,
        1.974,
        0.074,
    };

    static const std::array<TLetter, 256> kCharToLetter = {
        /* { 0,   */   0 /* } */,
        /* { 1,   */   0 /* } */,
        /* { 2,   */   0 /* } */,
        /* { 3,   */   0 /* } */,
        /* { 4,   */   0 /* } */,
        /* { 5,   */   0 /* } */,
        /* { 6,   */   0 /* } */,
        /* { 7,   */   0 /* } */,
        /* { 8,   */   0 /* } */,
        /* { 9,   */   0 /* } */,
        /* { 10,  */   0 /* } */,
        /* { 11,  */   0 /* } */,
        /* { 12,  */   0 /* } */,
        /* { 13,  */   0 /* } */,
        /* { 14,  */   0 /* } */,
        /* { 15,  */   0 /* } */,
        /* { 16,  */   0 /* } */,
        /* { 17,  */   0 /* } */,
        /* { 18,  */   0 /* } */,
        /* { 19,  */   0 /* } */,
        /* { 20,  */   0 /* } */,
        /* { 21,  */   0 /* } */,
        /* { 22,  */   0 /* } */,
        /* { 23,  */   0 /* } */,
        /* { 24,  */   0 /* } */,
        /* { 25,  */   0 /* } */,
        /* { 26,  */   0 /* } */,
        /* { 27,  */   0 /* } */,
        /* { 28,  */   0 /* } */,
        /* { 29,  */   0 /* } */,
        /* { 30,  */   0 /* } */,
        /* { 31,  */   0 /* } */,
        /* { 32,  */   0 /* } */,
        /* { 33,  */   0 /* } */,
        /* { 34,  */   0 /* } */,
        /* { 35,  */   0 /* } */,
        /* { 36,  */   0 /* } */,
        /* { 37,  */   0 /* } */,
        /* { 38,  */   0 /* } */,
        /* { 39,  */   0 /* } */,
        /* { 40,  */   0 /* } */,
        /* { 41,  */   0 /* } */,
        /* { 42,  */   0 /* } */,
        /* { 43,  */   0 /* } */,
        /* { 44,  */   0 /* } */,
        /* { 45,  */   0 /* } */,
        /* { 46,  */   0 /* } */,
        /* { 47,  */   0 /* } */,
        /* { 48,  */   0 /* } */,
        /* { 49,  */   0 /* } */,
        /* { 50,  */   0 /* } */,
        /* { 51,  */   0 /* } */,
        /* { 52,  */   0 /* } */,
        /* { 53,  */   0 /* } */,
        /* { 54,  */   0 /* } */,
        /* { 55,  */   0 /* } */,
        /* { 56,  */   0 /* } */,
        /* { 57,  */   0 /* } */,
        /* { 58,  */   0 /* } */,
        /* { 59,  */   0 /* } */,
        /* { 60,  */   0 /* } */,
        /* { 61,  */   0 /* } */,
        /* { 62,  */   0 /* } */,
        /* { 63,  */   0 /* } */,
        /* { 64,  */   0 /* } */,
        /* { 65,  */   1 /* } */,
        /* { 66,  */   2 /* } */,
        /* { 67,  */   3 /* } */,
        /* { 68,  */   4 /* } */,
        /* { 69,  */   5 /* } */,
        /* { 70,  */   6 /* } */,
        /* { 71,  */   7 /* } */,
        /* { 72,  */   8 /* } */,
        /* { 73,  */   9 /* } */,
        /* { 74,  */  10 /* } */,
        /* { 75,  */  11 /* } */,
        /* { 76,  */  12 /* } */,
        /* { 77,  */  13 /* } */,
        /* { 78,  */  14 /* } */,
        /* { 79,  */  15 /* } */,
        /* { 80,  */  16 /* } */,
        /* { 81,  */  17 /* } */,
        /* { 82,  */  18 /* } */,
        /* { 83,  */  19 /* } */,
        /* { 84,  */  20 /* } */,
        /* { 85,  */  21 /* } */,
        /* { 86,  */  22 /* } */,
        /* { 87,  */  23 /* } */,
        /* { 88,  */  24 /* } */,
        /* { 89,  */  25 /* } */,
        /* { 90,  */  26 /* } */,
        /* { 91,  */   0 /* } */,
        /* { 92,  */   0 /* } */,
        /* { 93,  */   0 /* } */,
        /* { 94,  */   0 /* } */,
        /* { 95,  */   0 /* } */,
        /* { 96,  */   0 /* } */,
        /* { 97,  */   1 /* } */,
        /* { 98,  */   2 /* } */,
        /* { 99,  */   3 /* } */,
        /* { 100, */   4 /* } */,
        /* { 101, */   5 /* } */,
        /* { 102, */   6 /* } */,
        /* { 103, */   7 /* } */,
        /* { 104, */   8 /* } */,
        /* { 105, */   9 /* } */,
        /* { 106, */  10 /* } */,
        /* { 107, */  11 /* } */,
        /* { 108, */  12 /* } */,
        /* { 109, */  13 /* } */,
        /* { 110, */  14 /* } */,
        /* { 111, */  15 /* } */,
        /* { 112, */  16 /* } */,
        /* { 113, */  17 /* } */,
        /* { 114, */  18 /* } */,
        /* { 115, */  19 /* } */,
        /* { 116, */  20 /* } */,
        /* { 117, */  21 /* } */,
        /* { 118, */  22 /* } */,
        /* { 119, */  23 /* } */,
        /* { 120, */  24 /* } */,
        /* { 121, */  25 /* } */,
        /* { 122, */  26 /* } */,
        /* { 123, */   0 /* } */,
        /* { 124, */   0 /* } */,
        /* { 125, */   0 /* } */,
        /* { 126, */   0 /* } */,
        /* { 127, */   0 /* } */,
        /* { 128, */   0 /* } */,
        /* { 129, */   0 /* } */,
        /* { 130, */   0 /* } */,
        /* { 131, */   0 /* } */,
        /* { 132, */   0 /* } */,
        /* { 133, */   0 /* } */,
        /* { 134, */   0 /* } */,
        /* { 135, */   0 /* } */,
        /* { 136, */   0 /* } */,
        /* { 137, */   0 /* } */,
        /* { 138, */   0 /* } */,
        /* { 139, */   0 /* } */,
        /* { 140, */   0 /* } */,
        /* { 141, */   0 /* } */,
        /* { 142, */   0 /* } */,
        /* { 143, */   0 /* } */,
        /* { 144, */   0 /* } */,
        /* { 145, */   0 /* } */,
        /* { 146, */   0 /* } */,
        /* { 147, */   0 /* } */,
        /* { 148, */   0 /* } */,
        /* { 149, */   0 /* } */,
        /* { 150, */   0 /* } */,
        /* { 151, */   0 /* } */,
        /* { 152, */   0 /* } */,
        /* { 153, */   0 /* } */,
        /* { 154, */   0 /* } */,
        /* { 155, */   0 /* } */,
        /* { 156, */   0 /* } */,
        /* { 157, */   0 /* } */,
        /* { 158, */   0 /* } */,
        /* { 159, */   0 /* } */,
        /* { 160, */   0 /* } */,
        /* { 161, */   0 /* } */,
        /* { 162, */   0 /* } */,
        /* { 163, */   0 /* } */,
        /* { 164, */   0 /* } */,
        /* { 165, */   0 /* } */,
        /* { 166, */   0 /* } */,
        /* { 167, */   0 /* } */,
        /* { 168, */   0 /* } */,
        /* { 169, */   0 /* } */,
        /* { 170, */   0 /* } */,
        /* { 171, */   0 /* } */,
        /* { 172, */   0 /* } */,
        /* { 173, */   0 /* } */,
        /* { 174, */   0 /* } */,
        /* { 175, */   0 /* } */,
        /* { 176, */   0 /* } */,
        /* { 177, */   0 /* } */,
        /* { 178, */   0 /* } */,
        /* { 179, */   0 /* } */,
        /* { 180, */   0 /* } */,
        /* { 181, */   0 /* } */,
        /* { 182, */   0 /* } */,
        /* { 183, */   0 /* } */,
        /* { 184, */   0 /* } */,
        /* { 185, */   0 /* } */,
        /* { 186, */   0 /* } */,
        /* { 187, */   0 /* } */,
        /* { 188, */   0 /* } */,
        /* { 189, */   0 /* } */,
        /* { 190, */   0 /* } */,
        /* { 191, */   0 /* } */,
        /* { 192, */   0 /* } */,
        /* { 193, */   0 /* } */,
        /* { 194, */   0 /* } */,
        /* { 195, */   0 /* } */,
        /* { 196, */   0 /* } */,
        /* { 197, */   0 /* } */,
        /* { 198, */   0 /* } */,
        /* { 199, */   0 /* } */,
        /* { 200, */   0 /* } */,
        /* { 201, */   0 /* } */,
        /* { 202, */   0 /* } */,
        /* { 203, */   0 /* } */,
        /* { 204, */   0 /* } */,
        /* { 205, */   0 /* } */,
        /* { 206, */   0 /* } */,
        /* { 207, */   0 /* } */,
        /* { 208, */   0 /* } */,
        /* { 209, */   0 /* } */,
        /* { 210, */   0 /* } */,
        /* { 211, */   0 /* } */,
        /* { 212, */   0 /* } */,
        /* { 213, */   0 /* } */,
        /* { 214, */   0 /* } */,
        /* { 215, */   0 /* } */,
        /* { 216, */   0 /* } */,
        /* { 217, */   0 /* } */,
        /* { 218, */   0 /* } */,
        /* { 219, */   0 /* } */,
        /* { 220, */   0 /* } */,
        /* { 221, */   0 /* } */,
        /* { 222, */   0 /* } */,
        /* { 223, */   0 /* } */,
        /* { 224, */   0 /* } */,
        /* { 225, */   0 /* } */,
        /* { 226, */   0 /* } */,
        /* { 227, */   0 /* } */,
        /* { 228, */   0 /* } */,
        /* { 229, */   0 /* } */,
        /* { 230, */   0 /* } */,
        /* { 231, */   0 /* } */,
        /* { 232, */   0 /* } */,
        /* { 233, */   0 /* } */,
        /* { 234, */   0 /* } */,
        /* { 235, */   0 /* } */,
        /* { 236, */   0 /* } */,
        /* { 237, */   0 /* } */,
        /* { 238, */   0 /* } */,
        /* { 239, */   0 /* } */,
        /* { 240, */   0 /* } */,
        /* { 241, */   0 /* } */,
        /* { 242, */   0 /* } */,
        /* { 243, */   0 /* } */,
        /* { 244, */   0 /* } */,
        /* { 245, */   0 /* } */,
        /* { 246, */   0 /* } */,
        /* { 247, */   0 /* } */,
        /* { 248, */   0 /* } */,
        /* { 249, */   0 /* } */,
        /* { 250, */   0 /* } */,
        /* { 251, */   0 /* } */,
        /* { 252, */   0 /* } */,
        /* { 253, */   0 /* } */,
        /* { 254, */   0 /* } */,
        /* { 255, */   0 /* } */,
    };

    TCode calcCode(const char * data, int n) {
        TCode res = 0;
        do { res <<= 5; res += kCharToLetter[*data++]; } while (--n > 0);
        return res;
    }

    bool loadFreqMap(const char * fname, TFreqMap & res) {
        auto & len = res.len;
        auto & prob = res.prob;

        len = 0;
        prob.clear();

        printf("[+] Loading n-gram file '%s'\n", fname);
        std::ifstream fin(fname);

        if (fin.good() == false) {
            printf("    Failed to open file\n");
            return false;
        }

        std::string gram;
        int32_t nfreq = 0;
        int64_t nTotal = 0;

        while (true) {
            fin >> gram >> nfreq;
            if (fin.eof()) break;

            nTotal += nfreq;
        }

        fin.clear();
        fin.seekg(0);

        while (true) {
            fin >> gram >> nfreq;
            if (fin.eof()) break;

            if (len == 0) {
                len = gram.size();
                prob.resize(1 << (5*len), 1.0);
            } else if (len != gram.size()) {
                printf("Error: loaded n-grams with vaying lengths\n");
                return false;
            }

            TCode idx = calcCode(gram.data(), len);
            if (prob[idx] < 0.5) {
                printf("Error: duplicate n-gram '%s'\n", gram.c_str());
                return false;
            }
            prob[idx] = std::log10(((double)(nfreq))/nTotal);
        }
        printf("    Total n-grams loaded = %g\n", (double) nTotal);

        res.pmin = std::log(100.01/nTotal);
        printf("    P-min = %g\n", res.pmin);
        for (auto & p : prob) {
            if (p >= 0.5) p = res.pmin;
        }

        return true;
    }

    bool generateClusters(const TParameters & params, int N, TClusters & clusters, const std::vector<int> & hint) {
        clusters.clear();
        clusters.resize(N);

        std::map<int, bool> used;
        for (int i = 0; i < N; ++i) {
            if (hint[i] < 0) continue;
            clusters[i] = hint[i];
            used[hint[i]] = true;
        }

        int cid = 0;
        for (int i = 0; i < N; ++i) {
            if (hint[i] >= 0) continue;
            //clusters[i] = rand()%params.maxClusters;
            while (used[cid]) {
                if (++cid >= params.maxClusters) cid = 0;
            }
            clusters[i] = cid;
            if (++cid >= params.maxClusters) cid = 0;
        }

        return true;
    }

    bool printClusterGoodness(const std::string & text, const TClusters & clusters) {
        int len = text.length();
        int ncorrect = 0;
        int ntotal = 0;
        for (int i = 0; i < len; ++i) {
            for (int j = i + 1; j < len; ++j) {
                bool ism0 = (text[i] == text[j]);
                bool ism1 = (clusters[i] == clusters[j]);
                if ((ism0 && ism1) || (!ism0 && !ism1)) {
                    ++ncorrect;
                }
                ++ntotal;
            }
        }

        printf("[+] Cluster goodness = %g\n", ((float)(ncorrect))/ntotal);

        return true;
    }

    float costF(const TSimilarityMap & ccMap, const TClusters & clusters) {
        int len = clusters.size();
        float res = 0.0f;
        for (int i = 0; i < len; ++i) {
            for (int j = i + 1; j < len; ++j) {
                if (clusters[i] == clusters[j]) {
                    res += 1.0f - ccMap.at(i).at(j);
                }
            }
        }

        return -res;
    }

    float costFUpdate(const TSimilarityMap & ccMap, const TClusters & clusters, int i, int cid, float c0) {
        int len = clusters.size();
        float res = -c0;
        for (int j = 0; j < len; ++j) {
            if (cid == clusters[j]) {
                res -= 1.0f - ccMap.at(i).at(j);
            }
            if (clusters[i] == clusters[j]) {
                res += 1.0f - ccMap.at(i).at(j);
            }
        }

        return -res;
    }

    TProb calcScore0(const TParameters & params, const TFreqMap & freqMap, const TClusters & txt, const TClusterToLetterMap & clMap) {
        TProb res = 0.0;

        auto n = txt.size();
        const auto & len  = freqMap.len;
        const auto & pmin = freqMap.pmin;
        const auto & prob = freqMap.prob;

        int nlet = 0;
        std::array<int, 26> letCount;
        letCount.fill(0);
        std::vector<TLetter> plain(n);
        for (int i = 0; i < n; ++i) {
            plain[i] = clMap.at(txt[i]);
            if (plain[i] > 0 && plain[i] <= 26) {
                ++letCount[plain[i] - 1];
                ++nlet;
            }
        }

        float letFreqCost = 0.0;
        for (int i = 0; i < 26; ++i) {
            float curf = 0.01*kEnglishLetterFreq[i] - ((float)(letCount[i]))/((float)(nlet));
            letFreqCost += curf*curf;
        }
        letFreqCost /= 26.0;
        letFreqCost = sqrt(letFreqCost);
        //printf("letFreqCost = %g\n", letFreqCost);

        int i1 = 0;
        int k = len;
        TCode curc = 0;
        TCode mask = (1 << 5*(len-1)) - 1;

        while (k > 0) {
            if (i1 >= n) return -1e100;
            auto c = plain[i1];
            ++i1;
            if (c > 0 && c <= 26) {
                curc <<= 5;
                curc += c;
                --k;
            } else {
                res += pmin;
            }
        }

        res += prob[curc];
        while (true) {
            curc &= mask;

            while (true) {
                if (i1 >= n) return res/n - params.wEnglishFreq*letFreqCost;
                auto c = plain[i1];
                ++i1;
                if (c > 0 && c <= 26) {
                    curc <<= 5;
                    curc += c;
                    break;
                } else {
                    res += pmin;
                }
            }

            res += prob[curc];
        }

        return res/n - params.wEnglishFreq*letFreqCost;
    }

    bool doSimulatedAnnealing3(
        const TParameters & params,
        const TSimilarityMap & ccMap,
        TClusters & clusters,
        const std::vector<int> & hint
        ) {

        int N = clusters.size();

        int ncc = 0;
        float ccavg = 0.0;
        for (int j = 0; j < N; ++j) {
            for (int i = j + 1; i < N; ++i) {
                ccavg += ccMap.at(j).at(i);
                ++ncc;
            }
        }
        ccavg /= ncc;
        printf("Average cc = %g\n", ccavg);

        float temp = params.temp0;
        float cost0 = costF(ccMap, clusters);

        auto clustersNew = clusters;
        auto clustersBest = clusters;

        bool hintMerge = false;

        int nRepeat = 0;
        float costOld = cost0;

        for (int iter = 0; iter < params.saMaxIterations; ++iter) {
            if (iter%5000 == 0) {
                printf("Iter %5d : temp = %16.4f, cost0 = %8.4f\n", iter, temp, cost0);
            }

            if (costOld == cost0) {
                if (++nRepeat == 0.1*params.saMaxIterations) {
                    printf("No improvement for long time. Aborting\n");
                    iter = 0.95*params.saMaxIterations + 1;
                }
            } else {
                nRepeat = 0;
                costOld = cost0;
            }

            bool doMerge = iter > 0.95*params.saMaxIterations;

            if (doMerge) {
                if (hintMerge == false) {
                    bool isValid = true;
                    std::map<TCluster, int> ch;
                    for (int k = 0; k < N; ++k) {
                        if (hint[k] < 0) continue;
                        int cid = clusters[k];
                        if (ch.find(cid) != ch.end()) {
                            if (ch[cid] != hint[k]) {
                                isValid = false;
                                break;
                            }
                        } else {
                            ch[cid] = hint[k];
                        }
                    }
                    printf("IsValid = %d\n", isValid);

                    for (int k0 = 0; k0 < N; ++k0) {
                        if (hint[k0] < 0) continue;
                        for (int k1 = k0 + 1; k1 < N; ++k1) {
                            if (hint[k1] < 0) continue;
                            if (hint[k0] != hint[k1]) continue;

                            if (clusters[k0] != clusters[k1]) {
                                printf("Merging clusters %d and %d based on hint\n", clusters[k0], clusters[k1]);

                                int cid0 = clusters[k0];
                                int cid1 = clusters[k1];

                                for (int i = 0; i < N; ++i) {
                                    if (clusters[i] == cid1) {
                                        clusters[i] = cid0;
                                    }
                                }
                            }
                        }
                    }
                    cost0 = costF(ccMap, clusters);
                    hintMerge = true;
                }

                clustersNew = clusters;

                int i0 = rand()%N;
                int i1 = rand()%N;
                while (clusters[i0] == clusters[i1] || ccMap.at(i0).at(i1) < ccavg) {
                    i0 = rand()%N;
                    i1 = rand()%N;
                }

                int cid0 = clusters[i0];
                int cid1 = clusters[i1];

                std::map<int, int> csize;
                for (const auto & c : clusters) {
                    ++csize[c];
                }
                int nonzero = csize.size();

                if (nonzero <= params.minClusters) {
                    printf("Reached min number of clusters. Done\n");
                    break;
                }

                break;

                int n0 = csize[cid0];
                int n1 = csize[cid1];

                --nonzero;

                for (int i = 0; i < N; ++i) {
                    if (clustersNew[i] == cid1) {
                        clustersNew[i] = cid0;
                    }
                }

                float cost1 = costF(ccMap, clustersNew);
                cost1 += 1.0f*n0*n1*(1.0f - ccavg);

                float delta = cost1 - cost0;
                if (delta > 0 || (std::exp(delta/temp) > frand())) {
                    printf("merge %d size %d      -  %d size %d\n", cid0, n0, cid1, n1);
                    cost0 = cost1 - 1.0f*n0*n1*(1.0f - ccavg);
                    clusters = clustersNew;
                }
            } else {
                clustersNew = clusters;
                clustersBest = clusters;

                float costBest = -1e10;
                float costCur = cost0;
                for (int k = 0; k < params.nChangePerIteration; ++k) {
                    int i = rand()%N;
                    while (hint[i] >= 0) {
                        i = rand()%N;
                    }

                    int cid = clustersNew[i];
                    while (cid == clustersNew[i]) {
                        clustersNew[i] = rand()%params.maxClusters;
                    }

                    costCur = costFUpdate(ccMap, clustersNew, i, cid, costCur);

                    if (costCur > costBest) {
                        costBest = costCur;
                        clustersBest = clustersNew;
                    }
                }

                float cost1 = costBest;

                float delta = cost1 - cost0;
                if (delta > 0 || (std::exp(delta/temp) > frand())) {
                    cost0 = cost1;
                    clusters = clustersBest;
                }
            }

            temp *= params.coolingRate;
        }

        return true;
    }

    bool doSimulatedAnnealing4(
        const TParameters & params,
        const TFreqMap & freqMap,
        const TClusters & clusters,
        TClusterToLetterMap & clMap,
        const std::vector<int> & hint
        ) {

        int N = clusters.size();

        clMap.clear();
        std::map<TLetter, bool> used;
        std::map<TCluster, bool> fixed;
        for (int i = 0; i < (int) hint.size(); ++i) {
            if (hint[i] < 0) continue;
            //if (used[hint[i]]) continue;
            if (fixed[clusters[i]]) continue;

            clMap[clusters[i]] = hint[i];
            used[hint[i]] = true;
            fixed[clusters[i]] = true;
            printf("Fixed %d\n", clusters[i]);
        }

        {
            int ii = 0;
            for (int i = 0; i < params.maxClusters; ++i) {
                if (fixed[i]) continue;
                while (used[ii]) ++ii;
                clMap[i] = ii++;
            }
        }

        float temp = params.temp0;
        float cost0 = calcScore0(params, freqMap, clusters, clMap);

        auto clMapSave = clMap;
        auto clMapBest = clMap;

        for (int iter = 0; iter < params.saMaxIterations; ++iter) {
            if (iter%5000 == 0) {
                printf("Iter %5d : temp = %16.4f, cost0 = %8.4f\n", iter, temp, cost0);
                for (int i = 0; i < N; ++i) {
                    if (clMap[clusters[i]] > 0 && clMap[clusters[i]] <= 26) {
                        printf("%c", 'a'+clMap[clusters[i]]-1);
                    } else {
                        printf(".");
                    }
                }
                printf("\n");
            }

            clMapSave = clMap;
            clMapBest = clMap;

            float costBest = -1e10;

            for (int k = 0; k < 1; ++k) {
                {
                    int i2 = rand()%params.maxClusters;
                    int i3 = rand()%params.maxClusters;
                    while (i2 == i3 || fixed[i2] || fixed[i3]) {
                        i2 = rand()%params.maxClusters;
                        i3 = rand()%params.maxClusters;
                    }
                    std::swap(clMap[i2], clMap[i3]);
                }

                int i2 = rand()%params.maxClusters;
                while (fixed[i2]) {
                    i2 = rand()%params.maxClusters;
                }
                int letterNew = rand()%27;
                while (clMap[i2] == letterNew) {
                    letterNew = rand()%27;
                }

                clMap[i2] = letterNew;

                float costCur = calcScore0(params, freqMap, clusters, clMap);
                if (costCur > costBest) {
                    costBest = costCur;
                    clMapBest = clMap;
                }
            }

            //float costBest = -1e10;
            //for (int k = 0; k < 1; ++k) {
            //    int i2 = rand()%params.maxClusters;
            //    int i3 = rand()%params.maxClusters;
            //    while (i2 == i3 || fixed[i2] || fixed[i3]) {
            //        i2 = rand()%params.maxClusters;
            //        i3 = rand()%params.maxClusters;
            //    }
            //    std::swap(clMap[i2], clMap[i3]);

            //    float costCur = calcScore0(params, freqMap, clusters, clMap);

            //    if (costCur > costBest) {
            //        costBest = costCur;
            //        clMapBest = clMap;
            //    }
            //}

            clMap = clMapBest;

            float cost1 = costBest;

            float delta = cost1 - cost0;
            if (delta > 0 || (std::exp(delta/temp) > frand())) {
                cost0 = cost1;
            } else {
                clMap = clMapSave;
            }

            temp *= params.coolingRate;
        }

        for (auto & cl : clMap) {
            printf("%d - %c\n", cl.first, 'a' + cl.second - 1);
        }

        for (int i = 0; i < 27; ++i) {
            printf("    %c - ", 'a' + i - 1);
            for (auto & cl : clMap) {
                if (cl.second == i) {
                    printf("%d ", cl.first);
                }
            }
            printf("\n");
        }

        return true;
    }

    bool doSimulatedAnnealing5(
        const TParameters & params,
        const TFreqMap & freqMap,
        const TSimilarityMap & ccMap,
        TClusters & clusters,
        TClusterToLetterMap & clMap,
        const std::vector<int> & hint
        ) {

        int N = clusters.size();

        int ncc = 0;
        float ccavg = 0.0;
        for (int j = 0; j < N; ++j) {
            for (int i = j + 1; i < N; ++i) {
                ccavg += ccMap.at(j).at(i);
                ++ncc;
            }
        }
        ccavg /= ncc;
        printf("Average cc = %g\n", ccavg);

        clMap.clear();
        std::map<TLetter, bool> used;
        std::map<TCluster, bool> fixed;
        for (int i = 0; i < (int) hint.size(); ++i) {
            if (hint[i] < 0) continue;
            //if (used[hint[i]]) continue;
            if (fixed[clusters[i]]) continue;

            clMap[clusters[i]] = hint[i];
            used[hint[i]] = true;
            fixed[clusters[i]] = true;
            printf("Fixed %d\n", clusters[i]);
        }

        {
            int ii = 0;
            for (int i = 0; i < params.maxClusters; ++i) {
                if (fixed[i]) continue;
                while (used[ii]) ++ii;
                clMap[i] = ii++;
            }
        }

        float wlm = params.wLanguageModel;
        float temp = params.temp0;
        float cost0CL = costF(ccMap, clusters);
        float cost0LM = calcScore0(params, freqMap, clusters, clMap);
        float cost0 = cost0CL + wlm*cost0LM;

        int nRepeat = 0;
        float costOld = cost0;

        auto tStart = std::chrono::high_resolution_clock::now();

        for (int iter = 0; iter < params.saMaxIterations; ++iter) {
            if (iter%10000 == 0) {
                auto tEnd = std::chrono::high_resolution_clock::now();
                float iterPerT = (float)(iter)/((float)(std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tStart).count()));

                printf("Speed: %6.2f iter/ms, Iter %5d : temp = %16.4f, cost0 = %8.4f, costCL = %8.4f, costLM = %8.4f\n",
                       iterPerT, iter, temp, cost0, cost0CL, cost0LM);
                for (int i = 0; i < clusters.size(); ++i) { printf("%c", 'a'+clMap[clusters[i]]-1); }
                printf("\n");
            }

            if (costOld == cost0) {
                //if (++nRepeat == 0.005*params.saMaxIterations) {
                //    printf("Heating\n");
                //    temp *= 100.1f;
                //}
                if (++nRepeat == 0.01*params.saMaxIterations) {
                    printf("No improvement for long time. Aborting\n");
                    iter = params.saMaxIterations + 1;
                }
            } else {
                nRepeat = 0;
                costOld = cost0;
            }

            int cid = -1;
            int i1l = -1;
            int i = rand()%N;
            int i1 = rand()%params.maxClusters;
            int i2 = rand()%params.maxClusters;
            int i3 = rand()%params.maxClusters;

            float costCurCL = cost0CL;
            float costCurLM = cost0LM;

            {
                while (hint[i] >= 0) {
                    i = rand()%N;
                }

                cid = clusters[i];
                while (cid == clusters[i]) {
                    clusters[i] = rand()%params.maxClusters;
                }

                costCurCL = costFUpdate(ccMap, clusters, i, cid, cost0CL);

                float cost1 = costCurCL + wlm*cost0LM;
                float delta = cost1 - cost0;
                if (delta > 0 || (std::exp(0.01f*(delta/temp)) > frand())) {
                    cost0CL = costCurCL;
                    cost0 = cost1;
                } else {
                    clusters[i] = cid;
                }
            }

            {
                while (i2 == i3 || fixed[i2] || fixed[i3]) {
                    i2 = rand()%params.maxClusters;
                    i3 = rand()%params.maxClusters;
                }
                std::swap(clMap[i2], clMap[i3]);

                while (fixed[i1]) {
                    i1 = rand()%params.maxClusters;
                }
                i1l = clMap[i1];
                int letterNew = rand()%27;
                while (clMap[i1] == letterNew) {
                    letterNew = rand()%27;
                }

                clMap[i1] = letterNew;

                costCurLM = calcScore0(params, freqMap, clusters, clMap);

                float cost1 = cost0CL + wlm*costCurLM;
                float delta = cost1 - cost0;
                if (delta > 0 || (std::exp(0.01f*(delta/temp)) > frand())) {
                    cost0LM = costCurLM;
                    cost0 = cost1;
                } else {
                    clMap[i1] = i1l;
                    std::swap(clMap[i2], clMap[i3]);
                }
            }

            temp *= params.coolingRate;
        }

        return true;
    }

    void getRandomCLMap(
        const TParameters & params,
        const TClusters & clusters,
        TClusterToLetterMap & clMap,
        const std::vector<int> & hint) {

        clMap.clear();
        std::map<TLetter, bool> used;
        std::map<TCluster, bool> fixed;
        for (int i = 0; i < (int) hint.size(); ++i) {
            if (hint[i] < 0) continue;
            //if (used[hint[i]]) continue;
            if (fixed[clusters[i]]) continue;

            clMap[clusters[i]] = hint[i];
            used[hint[i]] = true;
            fixed[clusters[i]] = true;
            printf("Fixed %d\n", clusters[i]);
        }

        {
            for (int i = 0; i < params.maxClusters; ++i) {
                if (fixed[i]) continue;
                clMap[i] = rand()%27;
            }
        }
    }

    bool subbreak(
        const TParameters & params,
        const TFreqMap & freqMap,
        const TClusters & clusters,
        TClusterToLetterMap & clMap,
        const std::vector<int> & hint) {
        int N = clusters.size();

        auto besta = clMap;
        getRandomCLMap(params, clusters, besta, hint);

        std::map<TLetter, bool> fixed;
        for (int i = 0; i < (int) hint.size(); ++i) {
            if (hint[i] < 0) continue;
            fixed[clusters[i]] = true;
        }

        float bestp = calcScore0(params, freqMap, clusters, besta);
        printf("Initial prob: %g\n", bestp);

        auto bestbestp = bestp;

        int nIters = 0;
        int nMainIters = params.nSubbreakIterations;
        while (nMainIters--) {
            if (++nIters > 10000) {
                getRandomCLMap(params, clusters, besta, hint);
                bestp = calcScore0(params, freqMap, clusters, besta);
                printf("reset\n");
                nIters = 0;
            }

            auto itera = besta;

            int nswaps = 3;
            for (int i = 0; i < nswaps; ++i) {
                int a0 = rand()%params.maxClusters;
                int a1 = rand()%params.maxClusters;
                while (a0 == a1) {
                    a0 = rand()%params.maxClusters;
                    a1 = rand()%params.maxClusters;
                }

                if (fixed[a0] || fixed[a1]) {
                } else {
                    std::swap(itera[a0], itera[a1]);
                }
            }

            auto iterp = calcScore0(params, freqMap, clusters, itera);
            auto cura = itera;
            for (int i = 0; i < 100; ++i) {
                int a0 = rand()%params.maxClusters;
                while (fixed[a0]) {
                    a0 = rand()%params.maxClusters;
                }

                int olda = cura[a0];
                cura[a0] = rand()%27;

                auto curp = calcScore0(params, freqMap, clusters, cura);
                if (curp > iterp) {
                    iterp = curp;
                    itera = cura;
                    i = 0;
                } else {
                    cura[a0] = olda;
                }
            }

            if (iterp > bestp) {
                besta = itera;
                bestp = iterp;

                if (bestp > bestbestp) {
                    bestbestp = bestp;
                    printf("[+] Best score = %g\n", bestp);
                    for (int i = 0; i < N; ++i) {
                        if (besta[clusters[i]] > 0 && besta[clusters[i]] <= 26) {
                            printf("%c", 'a'+besta[clusters[i]]-1);
                        } else {
                            printf(".");
                        }
                    }
                    printf("\n");
                }
                nIters = 0;
            }
        }

        clMap = besta;

        return true;
    }
}
