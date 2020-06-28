/*! \file subbreak2.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "subbreak2.h"

#include <array>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <cassert>
#include <algorithm>

namespace {

    static const std::array<float, 27> kEnglishLetterFreq = {
        0.000, // space
        8.167, // a
        1.492, // b
        2.782, // c
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

    static const std::array<float, 27> kEnglishLetterWithSpacesFreq = {
        11.505, // space
        7.228,  // a
        1.320,  // b
        2.462,  // c
        3.764,
        11.239,
        1.972,
        1.783,
        5.393,
        6.165,
        0.135,
        0.683,
        3.562,
        2.129,
        5.973,
        6.644,
        1.707,
        0.084,
        5.298,
        5.599,
        8.014,
        2.441,
        0.866,
        2.089,
        0.133,
        1.747,
        0.065,
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

// generate k random distinct ints in [0..n)
std::vector<int> subset(int k, int n) {
    std::vector<int> res(k);
    for (int i = 0; i < k; ++i) res[i] = i;

    for (int i = k; i < n; ++i) {
        const int j = rand()%(i + 1);
        if (j < k) {
            res[j] = i;
        }
    }

    return res;
}

template <typename T>
void shuffle(T & t, int start = -1, int end = -1, const Cipher::THint & hint = {}) {
    if (start == -1) start = 0;
    if (end == -1) end = t.size();

    for (int i = end - 1; i > start; --i) {
        int i0 = i;
        int i1 = rand()%(i - start + 1)+start;
        if (hint.size() > 0 && (hint[i0] != -1 || hint[i1] != -1)) {
            continue;
        }
        std::swap(t[i0], t[i1]);
    }
}

void calcHelpers(
        const Cipher::TParameters & params,
        const TClusterToLetterMap & clMap,
        const TClusters & clusters,
        Cipher::TPlainText & plainText,
        Cipher::TLetterCount & letterCount,
        Cipher::TClusterCount & clusterCount,
        Cipher::TClusterPos & clusterPos,
        int & nLetters,
        float & letterFreqCost) {

    int n = clusters.size();
    plainText.resize(n);
    letterCount.resize(n);
    std::fill(letterCount.begin(), letterCount.end(), 0);

    if (clusterCount.size() == 0) {
        clusterCount.resize(params.maxClusters, 0);
        clusterPos.resize(params.maxClusters);

        for (int i = 0; i < n; ++i) {
            ++clusterCount[clusters[i]];
            clusterPos[clusters[i]].push_back(i);
        }
    }

    nLetters = 0;
    for (int i = 0; i < n; ++i) {
        plainText[i] = clMap.at(clusters[i]);
        if (plainText[i] > 0 && plainText[i] <= 26) {
            ++letterCount[plainText[i]];
            ++nLetters;
        } else if (params.includeSpaces) {
            ++letterCount[0];
            ++nLetters;
        }
    }

    letterFreqCost = 0.0;
    {
        auto & freq = params.includeSpaces ? kEnglishLetterWithSpacesFreq : kEnglishLetterFreq;
        for (int i = 0; i < 27; ++i) {
            float curf = 0.01*freq[i] - ((float)(letterCount[i]))/((float)(nLetters));
            letterFreqCost += curf*curf;
        }
    }

    letterFreqCost /= 27.0;
    letterFreqCost = sqrt(letterFreqCost);
}

void updateHelpers(
        const Cipher::TParameters & params,
        const Cipher::TClusterCount & clusterCount,
        const Cipher::TClusterPos & clusterPos,
        const int nLetters,
        int c0, int c1,
        TClusterToLetterMap & clMap,
        Cipher::TPlainText & plainText,
        Cipher::TLetterCount & letterCount,
        float & letterFreqCost) {

    int p0 = clMap.at(c0);
    int p1 = clMap.at(c1);

    int nc0 = clusterCount[c0];
    int nc1 = clusterCount[c1];

    letterCount[p0] += nc1 - nc0;
    letterCount[p1] += nc0 - nc1;

    clMap.at(c0) = p1;
    clMap.at(c1) = p0;

    for (auto & i : clusterPos[c0]) plainText[i] = p1;
    for (auto & i : clusterPos[c1]) plainText[i] = p0;

    letterFreqCost = 0.0;
    {
        auto & freq = params.includeSpaces ? kEnglishLetterWithSpacesFreq : kEnglishLetterFreq;
        for (int i = 0; i < 27; ++i) {
            float curf = 0.01*freq[i] - ((float)(letterCount[i]))/((float)(nLetters));
            letterFreqCost += curf*curf;
        }
    }

    letterFreqCost /= 27.0;
    letterFreqCost = sqrt(letterFreqCost);
}

}

namespace Cipher {
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
        res.nTotal = 0;

        while (true) {
            fin >> gram >> nfreq;
            if (fin.eof()) break;

            res.nTotal += nfreq;
        }

        fin.clear();
        fin.seekg(0);

        while (true) {
            fin >> gram >> nfreq;
            if (fin.eof()) break;

            if (len == 0) {
                len = gram.size();
                prob.resize(1 << (5*len), 1.0);
            } else if (len != (int) gram.size()) {
                printf("Error: loaded n-grams with vaying lengths\n");
                return false;
            }

            TCode idx = calcCode(gram.data(), len);
            if (prob[idx] < 0.5) {
                printf("Error: duplicate n-gram '%s'\n", gram.c_str());
                return false;
            }
            prob[idx] = std::log10(((double)(nfreq))/res.nTotal);
        }
        printf("    Total n-grams loaded = %g\n", (double) res.nTotal);

        res.pmin = std::log(10000000.01/res.nTotal);
        printf("    P-min = %g\n", res.pmin);
        for (auto & p : prob) {
            if (p >= 0.5) p = res.pmin;
        }

        return true;
    }

    bool encryptExact(const TParameters & , const std::string & text, TClusters & clusters) {
        auto myCharToLetter = kCharToLetter;

		int k = 26;
        {
            for (auto & p : text) {
                if (myCharToLetter[p] == 0) {
                    myCharToLetter[p] = ++k;
                }
            }
            printf("[encryptExact] Unique symbols = %d\n", k);
        }

		std::vector<TClusterId> alphabetTransformation(k);
		for (int i = 0; i < (int) alphabetTransformation.size(); ++i) {
			alphabetTransformation[i] = i;
		}

		shuffle(alphabetTransformation);

		clusters.resize(text.size());

		for (int i = 0; i < (int) text.size(); ++i) {
			assert(myCharToLetter[text[i]] > 0);
			clusters[i] = alphabetTransformation[myCharToLetter[text[i]] - 1] + 1;
			printf("%d - %d\n", i, clusters[i]);
			assert(clusters[i] > 0 && clusters[i] <= k);
		}

        return true;
    }

	bool generateSimilarityMap(const TParameters & params, const std::string & text, TSimilarityMap & ccMap) {
		int n = text.size();

		ccMap.clear();
        ccMap.resize(n);
        for (auto & x : ccMap) x.resize(n);

		std::vector<float> waveformAccuracy(n);

		for (int i = 0; i < n; ++i) {
			float pError = frand();
			if (pError < params.waveformDetectionErrorP) {
				waveformAccuracy[i] = 1.0 - params.waveformDetectionErrorMin - frand()*(params.waveformDetectionErrorSig);
			} else {
				waveformAccuracy[i] = 1.0 - params.waveformDeviationMin - frand()*(params.waveformDeviationSig);
			}

			waveformAccuracy[i] += 2.0f*(0.5f - frand())*params.similarityNoiseSig;
			waveformAccuracy[i] = std::max(0.0f, waveformAccuracy[i]);
			waveformAccuracy[i] = std::min(1.0f, waveformAccuracy[i]);
		}

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				if (i == j) {
					ccMap[i][j].cc = 1.0;

					continue;
				}

				if (text[i] != text[j]) {
					float sim = params.similarityMismatchAvg + 2.0f*(0.5f - frand())*params.similarityMismatchSig;

					ccMap[i][j].cc = sim;
					ccMap[j][i].cc = sim;
				} else {
					float sim = waveformAccuracy[i]*waveformAccuracy[j];

					ccMap[i][j].cc = sim;
					ccMap[j][i].cc = sim;
				}
			}
		}

		return true;
	}

    bool generateClusters(const TParameters & params, int n, TClusters & clusters) {
        clusters.clear();
        clusters.resize(n);

        std::map<int, bool> used;
        for (int i = 0; i < n; ++i) {
            if (params.hint[i] < 0) continue;
            clusters[i] = params.hint[i];
            used[params.hint[i]] = true;
        }

        int cid = 0;
        for (int i = 0; i < n; ++i) {
            if (params.hint[i] >= 0) continue;
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
                    res += 1.0f - ccMap[i][j].cc;
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
                res -= 1.0f - ccMap[i][j].cc;
            }
            if (clusters[i] == clusters[j]) {
                res += 1.0f - ccMap[i][j].cc;
            }
        }

        return -res;
    }

    TProb calcScore0(const TParameters & params, const TFreqMap & freqMap, const TClusters & txt, const TClusterToLetterMap & clMap) {
        const int n = txt.size();
        const auto & len  = freqMap.len;
        const auto & prob = freqMap.prob;

        int nlet = 0;
        std::array<int, 27> letCount;
        letCount.fill(0);
        std::vector<TLetter> plain(n);
        for (int i = 0; i < n; ++i) {
            plain[i] = clMap.at(txt[i]);
            if (plain[i] > 0 && plain[i] <= 26) {
                ++letCount[plain[i]];
                ++nlet;
            } else if (params.includeSpaces) {
                ++letCount[0];
                ++nlet;
            }
        }

        float letFreqCost = 0.0;
        {
            auto & freq = params.includeSpaces ? kEnglishLetterWithSpacesFreq : kEnglishLetterFreq;
            for (int i = 0; i < 27; ++i) {
                float curf = 0.01*freq[i] - ((float)(letCount[i]))/((float)(nlet));
                letFreqCost += curf*curf;
            }
        }

        letFreqCost /= 27.0;
        letFreqCost = sqrt(letFreqCost);
        //printf("letFreqCost = %g\n", letFreqCost);

        int i1 = 0;
        int k = len;
        TCode curc = 0;
        TCode mask = (1 << 5*(len-1)) - 1;

        TProb res = 0.0;

        while (k > 0) {
            if (i1 >= n) return -1e100;
            auto c = plain[i1++];
            if (c > 0 && c <= 26) {
                curc <<= 5;
                curc += c;
                --k;
            } else {
                res += params.pNonAlphabetic;
            }
        }

        res += prob[curc];
        while (true) {
            curc &= mask;

            while (true) {
                if (i1 >= n) return res/n - params.wEnglishFreq*letFreqCost;
                auto c = plain[i1++];
                if (c > 0 && c <= 26) {
                    curc <<= 5;
                    curc += c;
                    break;
                } else {
                    res += params.pNonAlphabetic;
                }
            }

            res += prob[curc];
        }

        return res/n - params.wEnglishFreq*letFreqCost;
    }

    TProb calcScore1(
            const TParameters & params,
            const TFreqMap & freqMap,
            const std::vector<TLetter> & plain,
            const float letFreqCost) {
        TProb res = 0.0;

        const int n = plain.size();
        const auto & len  = freqMap.len;
        const auto & prob = freqMap.prob;

        //printf("letFreqCost = %g\n", letFreqCost);

        int i1 = 0;
        int k = len;
        TCode curc = 0;
        TCode mask = (1 << 5*(len-1)) - 1;

        while (k > 0) {
            if (i1 >= n) return -1e100;
            auto c = plain[i1++];
            if (c > 0 && c <= 26) {
                curc <<= 5;
                curc += c;
                --k;
            } else {
                res += params.pNonAlphabetic;
            }
        }

        res += prob[curc];
        while (true) {
            curc &= mask;

            while (true) {
                if (i1 >= n) return res/n - params.wEnglishFreq*letFreqCost;
                auto c = plain[i1++];
                if (c > 0 && c <= 26) {
                    curc <<= 5;
                    curc += c;
                    break;
                } else {
                    res += params.pNonAlphabetic;
                }
            }

            res += prob[curc];
        }

        return res/n - params.wEnglishFreq*letFreqCost;
    }

    bool doSimulatedAnnealing3(
        const TParameters & params,
        const TSimilarityMap & ccMap,
        TClusters & clusters
        ) {

        int n = clusters.size();

        int ncc = 0;
        float ccavg = 0.0;
        for (int j = 0; j < n; ++j) {
            for (int i = j + 1; i < n; ++i) {
                ccavg += ccMap[j][i].cc;
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
                    std::map<TClusterId, int> ch;
                    for (int k = 0; k < n; ++k) {
                        if (params.hint[k] < 0) continue;
                        int cid = clusters[k];
                        if (ch.find(cid) != ch.end()) {
                            if (ch[cid] != params.hint[k]) {
                                isValid = false;
                                break;
                            }
                        } else {
                            ch[cid] = params.hint[k];
                        }
                    }
                    printf("IsValid = %d\n", isValid);

                    for (int k0 = 0; k0 < n; ++k0) {
                        if (params.hint[k0] < 0) continue;
                        for (int k1 = k0 + 1; k1 < n; ++k1) {
                            if (params.hint[k1] < 0) continue;
                            if (params.hint[k0] != params.hint[k1]) continue;

                            if (clusters[k0] != clusters[k1]) {
                                printf("Merging clusters %d and %d based on hint\n", clusters[k0], clusters[k1]);

                                int cid0 = clusters[k0];
                                int cid1 = clusters[k1];

                                for (int i = 0; i < n; ++i) {
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

                int i0 = rand()%n;
                int i1 = rand()%n;
                while (clusters[i0] == clusters[i1] || ccMap[i0][i1].cc < ccavg) {
                    i0 = rand()%n;
                    i1 = rand()%n;
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

                for (int i = 0; i < n; ++i) {
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
                    int i = rand()%n;
                    while (params.hint[i] >= 0) {
                        i = rand()%n;
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
        TClusterToLetterMap & clMap
        ) {

        int n = clusters.size();

        clMap.clear();
        std::map<TLetter, bool> used;
        std::map<TClusterId, bool> fixed;
        for (int i = 0; i < (int) params.hint.size(); ++i) {
            if (params.hint[i] < 0) continue;
            //if (used[params.hint[i]]) continue;
            if (fixed[clusters[i]]) continue;

            clMap[clusters[i]] = params.hint[i];
            used[params.hint[i]] = true;
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
                for (int i = 0; i < n; ++i) {
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
        TClusterToLetterMap & clMap
        ) {

        int n = clusters.size();

        int ncc = 0;
        float ccavg = 0.0;
        for (int j = 0; j < n; ++j) {
            for (int i = j + 1; i < n; ++i) {
                ccavg += ccMap[j][i].cc;
                ++ncc;
            }
        }
        ccavg /= ncc;
        printf("Average cc = %g\n", ccavg);

        clMap.clear();
        std::map<TLetter, bool> used;
        std::map<TClusterId, bool> fixed;
        for (int i = 0; i < (int) params.hint.size(); ++i) {
            if (params.hint[i] < 0) continue;
            //if (used[params.hint[i]]) continue;
            if (fixed[clusters[i]]) continue;

            clMap[clusters[i]] = params.hint[i];
            used[params.hint[i]] = true;
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
                for (int i = 0; i < (int) clusters.size(); ++i) { printf("%c", 'a'+clMap[clusters[i]]-1); }
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
            int i = rand()%n;
            int i1 = rand()%params.maxClusters;
            int i2 = rand()%params.maxClusters;
            int i3 = rand()%params.maxClusters;

            float costCurCL = cost0CL;
            float costCurLM = cost0LM;

            {
                while (params.hint[i] >= 0) {
                    i = rand()%n;
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
        TClusterToLetterMap & clMap) {

        clMap.clear();

        for (int i = 0; i < params.maxClusters; ++i) {
            clMap[i] = rand()%27;
        }

        std::map<int, std::vector<int>> options;
        for (int i = 0; i < (int) params.hint.size(); ++i) {
            if (params.hint[i] < 0) continue;
            options[clusters[i]].push_back(params.hint[i]);
        }

        for (auto & option : options) {
            int n = option.second.size();
            clMap[option.first] = option.second[rand()%n];
        }
    }

    bool subbreak(
        const TParameters & params,
        const TFreqMap & freqMap,
        TResult & result) {
        const auto & clusters = result.clusters;
        auto & clMap = result.clMap;
        auto & minBestP = result.p;

        std::map<TLetter, bool> fixed;
        for (int i = 0; i < (int) params.hint.size(); ++i) {
            if (params.hint[i] < 0) continue;
            fixed[clusters[i]] = true;
        }

        auto besta = clMap;
        getRandomCLMap(params, clusters, besta);
        float bestp = calcScore0(params, freqMap, clusters, besta);

        int nIters = 0;
        int nMainIters = params.nSubbreakIterations;
        while (nMainIters--) {
            if (++nIters > 100) {
                getRandomCLMap(params, clusters, besta);
                bestp = calcScore0(params, freqMap, clusters, besta);
                nIters = 0;
            }

            auto itera = besta;

            int nswaps = 3;
            for (int i = 0; i < nswaps; ++i) {
                int a0 = rand()%params.maxClusters;
                int a1 = rand()%params.maxClusters;
                while (a0 == a1 || fixed[a0] || fixed[a1]) {
                    a0 = rand()%params.maxClusters;
                    a1 = rand()%params.maxClusters;
                }

                std::swap(itera[a0], itera[a1]);
            }

            auto iterp = calcScore0(params, freqMap, clusters, itera);
            auto cura = itera;
            for (int i = 0; i < 5; ++i) {
                int a0 = rand()%params.maxClusters;
                int a1 = rand()%params.maxClusters;
                while (a0 == a1 || fixed[a0] || fixed[a1]) {
                    a0 = rand()%params.maxClusters;
                    a1 = rand()%params.maxClusters;
                }

                std::swap(cura[a0], cura[a1]);

                auto curp = calcScore0(params, freqMap, clusters, cura);
                if (curp > iterp) {
                    iterp = curp;
                    itera = cura;
                    i = 0;
                } else {
					std::swap(cura[a0], cura[a1]);
                }
            }

            if (iterp > bestp) {
                besta = itera;
                bestp = iterp;

                if (bestp > minBestP) {
                    minBestP = bestp;
                    clMap = besta;
                }
                nIters = 0;
            }
        }

        return true;
    }

    bool subbreak1(
        const TParameters & params,
        const TFreqMap & freqMap,
        TResult & result) {

        // helpers
        TPlainText plainText;
        TLetterCount letterCount;
        TClusterCount clusterCount;
        TClusterPos clusterPos;

        int nLetters = 0;
        float letterFreqCost = 0.0f;

        const auto & clusters = result.clusters;
        auto & clMap = result.clMap;
        auto & minBestP = result.p;

        auto besta = clMap;
        getRandomCLMap(params, clusters, besta);

        std::map<TLetter, bool> fixed;
        for (int i = 0; i < (int) params.hint.size(); ++i) {
            if (params.hint[i] < 0) continue;
            fixed[clusters[i]] = true;
        }

        ::calcHelpers(params, besta, clusters, plainText, letterCount, clusterCount, clusterPos, nLetters, letterFreqCost);
        float bestp = calcScore1(params, freqMap, plainText, letterFreqCost);

        int nIters = 0;
        int nMainIters = params.nSubbreakIterations;
        while (nMainIters--) {
            if (++nIters > 100) {
                getRandomCLMap(params, clusters, besta);

                ::calcHelpers(params, besta, clusters, plainText, letterCount, clusterCount, clusterPos, nLetters, letterFreqCost);
                bestp = calcScore1(params, freqMap, plainText, letterFreqCost);

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
                    //std::swap(itera[a0], itera[a1]);
                    ::updateHelpers(params, clusterCount, clusterPos, nLetters, a0, a1, itera, plainText, letterCount, letterFreqCost);
                }
            }

            auto iterp = calcScore1(params, freqMap, plainText, letterFreqCost);

            auto cura = itera;
            for (int i = 0; i < 10; ++i) {
                int a0 = rand()%params.maxClusters;
                int a1 = rand()%params.maxClusters;
                while (a0 == a1) {
                    a0 = rand()%params.maxClusters;
                    a1 = rand()%params.maxClusters;
                }

                //std::swap(cura[a0], cura[a1]);
                ::updateHelpers(params, clusterCount, clusterPos, nLetters, a0, a1, cura, plainText, letterCount, letterFreqCost);
                auto curp = calcScore1(params, freqMap, plainText, letterFreqCost);

                if (curp > iterp) {
                    iterp = curp;
                    itera = cura;
                    i = 0;
                } else {
					//std::swap(cura[a0], cura[a1]);
                    ::updateHelpers(params, clusterCount, clusterPos, nLetters, a0, a1, cura, plainText, letterCount, letterFreqCost);
                }
            }

            if (iterp > bestp) {
                besta = itera;
                bestp = iterp;

                if (bestp > minBestP) {
                    minBestP = bestp;
                    clMap = besta;
                }
                nIters = 0;
            }
        }

        return true;
    }

    bool generateClustersInitialGuess(
            const TParameters & params,
            const TSimilarityMap & ccMap,
            TClusters & clusters) {
        int n = ccMap.size();

        int nClusters = n;
        clusters.resize(n);
        for (int i = 0; i < n; ++i) {
            clusters[i] = i;
        }

        struct Pair {
            int i;
            int j;
            double cc;

            bool operator < (const Pair & a) const { return cc > a.cc; }
        };

        std::vector<Pair> ccPairs;
        for (int i = 0; i < n - 1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                ccPairs.emplace_back(Pair{i, j, ccMap[i][j].cc});
            }
        }

        std::sort(ccPairs.begin(), ccPairs.end());

        {
            std::vector<bool> used(n);

            for (int k = 0; k < (int) ccPairs.size(); ++k) {
                int i = ccPairs[k].i;
                int j = ccPairs[k].j;

                if (clusters[i] == clusters[j]) continue;

                int cidi = clusters[i];
                int cidj = clusters[j];

                if (used[cidi] || used[cidj]) continue;

                if (cidi > cidj) {
                    std::swap(cidi, cidj);
                }

                for (int p = 0; p < n; ++p) {
                    if (clusters[p] == cidj) {
                        clusters[p] = cidi;
                    }
                }
                used[cidj] = true;
                --nClusters;

                if (nClusters <= params.maxClusters) break;
            }
        }

        {
            int cnt = 0;
            std::map<int, int> used;
            for (auto & cid : clusters) {
                if (used[cid] > 0) continue;
                used[cid] = ++cnt;
            }

            for (auto & cid : clusters) {
                cid = used[cid] - 1;
            }
        }

        printf("nClusters = %d\n", nClusters);
        for (auto & cid : clusters) {
            assert(cid >= 0 && cid < params.maxClusters);
        }

        return true;
    }

    bool mutateClusters(const TParameters & params, TClusters & clusters) {
        int n = clusters.size();

        //for (int i = 0; i < 3; ++i) {
        //    int i0 = rand()%n;
        //    int i1 = rand()%n;
        //    std::swap(clusters[i0], clusters[i1]);
        //}

        for (int i = 0; i < 1; ++i) {
            int idx = rand()%n;
            clusters[idx] = rand()%params.maxClusters;
            //auto p = clusters[idx];
            //while (clusters[idx] == p) {
            //    clusters[idx] = rand()%params.maxClusters;
            //}
        }

        //std::map<int, bool> used;
        //for (int i = 0; i < n; ++i) {
        //    if (hint[i] < 0) continue;
        //    clusters[i] = hint[i];
        //    used[hint[i]] = true;
        //}

        return true;
    }

    double calcPClusters(
            const TParameters & ,
            const TSimilarityMap & ,
            const TSimilarityMap & logMap,
            const TSimilarityMap & logMapInv,
            const TClusters & clusters,
            const TClusterToLetterMap &) {

        double res = 0.0;
        int n = clusters.size();

        for (int j = 0; j < n - 1; ++j) {
            for (int i = j + 1; i < n; ++i) {
                if (clusters[i] == clusters[j]) {
                //if (clMap.at(clusters[i]) == clMap.at(clusters[j])) {
                    res += logMap[j][i].cc;
                } else {
                    res += logMapInv[j][i].cc;
                }
            }
        }

        //if (clMap.empty()) {
        //    for (int j = 0; j < n - 1; ++j) {
        //        for (int i = j + 1; i < n; ++i) {
        //            if (clusters[i] == clusters[j]) {
        //                res += logMap[j][i].cc;
        //            } else {
        //                res += logMapInv[j][i].cc;
        //            }
        //        }
        //    }
        //    res *= 2.0;
        //} else {
        //    for (int j = 0; j < n - 1; ++j) {
        //        for (int i = j + 1; i < n; ++i) {
        //            if (clMap.at(clusters[i]) == clMap.at(clusters[j])) {
        //                res += logMap[j][i].cc;
        //            } else {
        //                res += logMapInv[j][i].cc;
        //            }
        //        }
        //    }
        //}

        return res;
    }

    bool normalizeSimilarityMap(
            const TParameters & ,
            TSimilarityMap & ccMap,
            TSimilarityMap & logMap,
            TSimilarityMap & logMapInv) {
        int n = ccMap.size();

        double ccMin = std::numeric_limits<double>::max();
        double ccMax = std::numeric_limits<double>::min();

        for (int j = 0; j < n - 1; ++j) {
            for (int i = j + 1; i < n; ++i) {
                ccMin = std::min(ccMin, ccMap[j][i].cc);
                ccMax = std::max(ccMax, ccMap[j][i].cc);
            }
        }

        ccMin -= 1e-6;
        ccMax += 1e-6;

        printf("ccMax = %g, ccMin = %g\n", ccMax, ccMin);

        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                if (i == j) {
                    ccMap[j][i].cc = 1.0;
                    continue;
                }

                auto & v = ccMap[j][i].cc;
                v = (v - ccMin)/(ccMax - ccMin);
                //v = 1.0 - std::exp(-1.1f*v);
            }
        }

        logMap = ccMap;

        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                if (i == j) {
                    logMap[j][i].cc = 0.0;
                    continue;
                }

                auto & v = logMap[j][i].cc;
                v = std::log(v);
            }
        }

        logMapInv = ccMap;

        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                if (i == j) {
                    logMap[j][i].cc = -1e6;
                    continue;
                }

                auto & v = logMapInv[j][i].cc;
                v = std::log(1.0 - v);
            }
        }

        return true;
    }

    void printText(const TClusters & t) {
        for (auto & c : t) {
            if (c >= 1 && c <= 26) printf("%c", 'a' + c - 1); else printf(".");
        }
        printf("\n");
    }

    void printText(const TClusters & t, const TClusterToLetterMap & clMap) {
        for (auto & cc : t) {
			auto c = clMap.at(cc);
            if (c >= 1 && c <= 26) printf("%c", 'a' + c - 1); else printf(".");
        }
        printf("\n");
    }

    //
    // Processor
    //

    Processor::Processor() {
    }

    bool Processor::init(
            const TParameters & params,
            const TFreqMap & freqMap,
            const TSimilarityMap & similarityMap) {
        m_params = params;
        m_freqMap = &freqMap;
        m_similarityMap = similarityMap;

        normalizeSimilarityMap(m_params, m_similarityMap, m_logMap, m_logMapInv);
        generateClustersInitialGuess(m_params, m_similarityMap, m_curResult.clusters);

        auto clustersNew = m_curResult.clusters;

        Cipher::mutateClusters(m_params, clustersNew);
        m_pCur = calcPClusters(m_params, m_similarityMap, m_logMap, m_logMapInv, clustersNew, m_curResult.clMap);

        return true;
    }

    bool Processor::setPNonAlphabetic(TProb p) {
        m_params.pNonAlphabetic = p;

        return true;
    }

    bool Processor::setWEnglishFreq(float w) {
        m_params.wEnglishFreq = w;
        m_curResult.p *= 1.001;
        m_pCur *= 1.001;

        return true;
    }

    bool Processor::setHint(const THint & hint) {
        m_params.hint = hint;

        return true;
    }

    bool Processor::compute() {
        auto clustersNew = m_curResult.clusters;

        for (int iter = 0; iter < m_params.nMHIters; ++iter) {
            clustersNew = m_curResult.clusters;
            Cipher::mutateClusters(m_params, clustersNew);
            auto pNew = calcPClusters(m_params, m_similarityMap, m_logMap, m_logMapInv, clustersNew, m_curResult.clMap);

            auto u = frand();
            //auto alpha = pNew/pCur;
            auto alpha = pNew > m_pCur ? 1.0 : std::exp((pNew - m_pCur));

            //printf("pNew = %g, pCur = %g, alpha = %g\n", pNew, pCur, alpha);

            if (u <= alpha) {
                m_curResult.clusters = clustersNew;
                m_pCur = pNew;

                if (++m_nMHInitialIters > m_params.nMHInitialIters) {
                    m_curResult.p *= 1.0001;
                    auto saveHint = m_params.hint;
                    for (auto & hint : m_params.hint) {
                        if (hint < 0) continue;
                        if (frand() > 0.10) {
                            hint = -1;
                        }
                    }
                    Cipher::subbreak(m_params, *m_freqMap, m_curResult);
                    m_params.hint = saveHint;
                }
            } else {
                m_pCur *= 1.00001;
                //auto x = m_params.nSubbreakIterations;
                //m_params.nSubbreakIterations = 100;
                //Cipher::subbreak(m_params, *m_freqMap, m_curResult);
                //m_params.nSubbreakIterations = x;
            }

            ++m_curResult.id;
        }

        return true;
    }

    const TResult & Processor::getResult() const {
        return m_curResult;
    }

    const TSimilarityMap & Processor::getSimilarityMap() const {
        return m_similarityMap;
    }

}
