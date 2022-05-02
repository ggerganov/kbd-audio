/*! \file subbreak3.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "subbreak3.h"
#include "constants.h"

#include <array>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <cassert>
#include <algorithm>

namespace {

    static const std::array<float, 28> kEnglishLetterWithSpacesFreq = {
        0.000,  // unused
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
        11.505, // space
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
        /* { 32,  */  27 /* } */,
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
        /* { 95,  */  27 /* } */,
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

}

namespace Cipher {
    TCode calcCode(const char * data, int n) {
        TCode res = 0;
        do { res <<= 5; res += kCharToLetter[*data++]; } while (--n > 0);
        return res;
    }

    bool loadFreqMap(const char * fname, TFreqMap & res, double pmin) {
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
        int64_t nfreq = 0;
        res.nTotal = 0;

        std::unordered_map<TCode, int64_t> pi64;

        while (true) {
            fin >> gram >> nfreq;
            if (fin.eof()) break;

            if (len == 0) {
                len = gram.size();
            } else if (len != (int) gram.size()) {
                printf("Error: loaded n-grams with vaying lengths\n");
                return false;
            }

            TCode idx = calcCode(gram.data(), len);
            if (pi64[idx] > 0) {
                printf("Error: duplicate n-gram '%s'\n", gram.c_str());
                return false;
            }
            pi64[idx] = nfreq;
            res.nTotal += nfreq;
        }
        printf("    Total n-grams loaded = %g\n", (double) res.nTotal);

        // compute wildcard frequencies
        {
            const auto tStart = std::chrono::steady_clock::now();

            std::vector<int32_t> masks((1 << len) - 1);
            for (int m = 0; m < (int) masks.size(); ++m) {
                auto & mask = masks[m];

                mask = 0;
                int32_t l = 0b11111;
                int32_t mm = m;
                for (int i = 0; i < len; ++i) {
                    if (mm & 1) {
                        mask |= l;
                    }
                    mm >>= 1;
                    l <<= 5;
                }

                // print mask in binary
                printf("    Mask %2d: ", m);
                for (int i = 0; i < 5*len; ++i) {
                    printf("%d", (mask >> (5*len - 1 - i)) & 1);
                }
                printf("\n");
            }

            {
                std::unordered_map<TCode, int64_t> pwild;
                for (const auto & [i, p] : pi64) {
                    if (p == 0) {
                        printf("i = %d, p == 0 - should not happen\n", i);
                        return false;
                    }

                    for (const auto & mask : masks) {
                        pwild[i & mask] += p;
                    }
                }

                printf("Size of pwild = %d\n", (int) pwild.size());

                for (auto & [i, p] : pwild) {
                    pi64[i] = p;
                }

                printf("Size of pi64 = %d\n", (int) pi64.size());
            }

            if (pi64[0] != res.nTotal) {
                printf("Error: wildcard probability mismatch - p[0] = %ld, expected %ld\n", (long) pi64[0], (long) res.nTotal);
                return false;
            }

            const auto tEnd = std::chrono::steady_clock::now();
            const auto tDiff = std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tStart).count();
            printf("    Wildcard probabilities computed in %g ms\n", (double) tDiff);
        }

        // compute probabilities
        {
            const auto tStart = std::chrono::steady_clock::now();

            res.pmin = std::log10(pmin);
            printf("    P-min = %g\n", res.pmin);

            for (auto & [i, p] : pi64) {
                if (p == 0) {
                    printf("i = %d, p == 0 - should not happen\n", i);
                    return false;
                } else {
                    double pp = double(p)/res.nTotal;
                    if (pp < pmin) {
                        prob[i] = res.pmin;
                    } else {
                        prob[i] = std::log10(pp);
                    }
                }
            }

            const auto tEnd = std::chrono::steady_clock::now();
            const auto tDiff = std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tStart).count();
            printf("    Probability computation time = %g ms\n", (double) tDiff);
        }

        return true;
    }

    bool saveFreqMapBinary(const char * fname, const TFreqMap & freqMap) {
        std::ofstream fout(fname, std::ios::binary);
        if (fout.good() == false) {
            printf("    Failed to open file '%s'\n", fname);
            return false;
        }

        fout.write((const char *) &freqMap.len,    sizeof(freqMap.len));
        fout.write((const char *) &freqMap.nTotal, sizeof(freqMap.nTotal));
        fout.write((const char *) &freqMap.pmin,   sizeof(freqMap.pmin));

        { int32_t n = freqMap.prob.size(); fout.write((const char *) &n, sizeof(n)); }

        {
            std::map<TCode, TProb> sorted(freqMap.prob.begin(), freqMap.prob.end());

            std::vector<TCode> keys;
            for (const auto & [i, p] : sorted) {
                if (p == freqMap.pmin) {
                    continue;
                }
                keys.push_back(i);
            }

            int32_t is = 0;
            while (is < (int) keys.size()) {
                fout.write((const char *) &keys[is],         sizeof(keys[is]));
                fout.write((const char *) &sorted[keys[is]], sizeof(sorted[keys[is]]));

                if (is == (int) keys.size() - 1) break;

                int32_t ie = is + 1;
                while (keys[ie] - keys[ie - 1] < 256) {
                    ++ie;
                    if (ie == (int) keys.size()) break;
                    if (ie == is + 255) break;
                }

                uint8_t n = (uint8_t) (ie - is - 1);
                fout.write((const char *) &n, sizeof(n));
                if (n > 0) {
                    for (int i = is + 1; i < ie; ++i) {
                        uint8_t d = (uint8_t) (keys[i] - keys[i - 1]);
                        fout.write((const char *) &d,               sizeof(d));
                        fout.write((const char *) &sorted[keys[i]], sizeof(sorted[keys[i]]));
                    }
                }

                is = ie;
            }
        }

        return true;
    }

    bool loadFreqMapBinary(const char * fname, TFreqMap & freqMap) {
        std::ifstream fin(fname, std::ios::binary);
        if (fin.good() == false) {
            printf("    Failed to open file '%s'\n", fname);
            return false;
        }

        fin.read((char *) &freqMap.len,    sizeof(freqMap.len));
        fin.read((char *) &freqMap.nTotal, sizeof(freqMap.nTotal));
        fin.read((char *) &freqMap.pmin,   sizeof(freqMap.pmin));

        {
            int32_t n;
            fin.read((char *) &n, sizeof(n));

            TCode curi;
            TProb curp;

            while (n > 0) {
                fin.read((char *) &curi, sizeof(curi));
                fin.read((char *) &curp, sizeof(curp));
                freqMap.prob[curi] = curp;

                uint8_t n8;
                fin.read((char *) &n8, sizeof(n8));

                if (n8 > 0) {
                    for (int i = 0; i < n8; ++i) {
                        uint8_t d;
                        fin.read((char *) &d,    sizeof(d));
                        fin.read((char *) &curp, sizeof(curp));
                        freqMap.prob[curi + d] = curp;
                        curi += d;
                    }

                    n -= n8;
                }

                --n;
            }
        }

        return true;
    }

    bool encryptExact(const TParameters & , const std::string & text, TClusters & clusters) {
        auto myCharToLetter = kCharToLetter;

		int k = 27;
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
			assert(clusters[i] > 0 && clusters[i] <= k);
		}

        return true;
    }

    TProb calcScore(
            const TParameters & params,
            const TFreqMap & freqMap,
            const std::vector<TLetter> & plain,
                  std::vector<TProb> & memo) {
        TProb res = 0.0;

        const int n = plain.size();
        const auto & len  = freqMap.len;
        const auto & prob = freqMap.prob;

        int nlet = 0;
        std::array<int, 28> letCount;
        letCount.fill(0);
        for (int i = 0; i < n; ++i) {
            if (plain[i] >= 0 && plain[i] <= 27) {
                ++letCount[plain[i]];
                ++nlet;
            }
        }

        float letFreqCost = 0.0;
        {
            auto & freq = kEnglishLetterWithSpacesFreq;
            for (int i = 0; i <= 27; ++i) {
                float curf = 0.01*freq[i] - ((float)(letCount[i]))/((float)(nlet));
                letFreqCost += curf*curf;
            }
        }

        letFreqCost /= 28.0;
        letFreqCost = sqrt(letFreqCost);

        int i1 = 0;
        int k = len;
        TCode curc = 0;
        TCode mask = (1 << 5*(len-1)) - 1;

        while (k > 0) {
            if (i1 >= n) return -1e100;
            auto c = plain[i1++];
            curc <<= 5;
            curc += c;
            --k;
        }

        if (memo[i1 - 1] > 0.5) {
            memo[i1 - 1] = prob.find(curc) == prob.end() ? freqMap.pmin : prob.at(curc);
        }
        res += memo[i1 - 1];

        while (true) {
            curc &= mask;

            auto c = plain[i1++];
            curc <<= 5;
            curc += c;

            if (memo[i1 - 1] > 0.5) {
                memo[i1 - 1] = prob.find(curc) == prob.end() ? freqMap.pmin : prob.at(curc);
            }
            res += memo[i1 - 1];

            if (i1 >= n) break;
        }

        return res/n - params.wEnglishFreq*letFreqCost;
    }

    TClusterToLetterMap getNullCLMap(const TClusters & clusters) {
        TClusterToLetterMap result;

        for (int i = 0; i < (int) clusters.size(); ++i) {
            result[clusters[i]] = 0;
        }

        return result;
    }

    void translate(
            const TClusterToLetterMap & clMap,
            const TClusters & clusters,
            std::vector<TLetter> & plain) {
        plain.clear();
        for (int i = 0; i < (int) clusters.size(); ++i) {
            plain.push_back(clMap.at(clusters[i]));
        }
    }

    bool beamSearch(
        const TParameters & params,
        const TFreqMap & freqMap,
        TResult & result) {
        const auto & clusters = result.clusters;

        struct THypothesis {
            TProb p;
            TClusterToLetterMap clMap;
            std::vector<TLetter> plain;
            std::vector<int> nused;
            std::vector<TProb> memo;
        };

        const int N = clusters.size();
        const int nSymbols = 27;
        const int nHypothesesToKeep = params.nHypothesesToKeep;

        int nCur = 0;
        std::vector<THypothesis> hypothesesCur;
        std::vector<THypothesis> hypothesesNew;
        hypothesesCur.resize(nHypothesesToKeep);
        hypothesesNew.resize(nHypothesesToKeep*nSymbols);

        {
            auto & hcur = hypothesesCur[0];
            hcur = { 0.0, getNullCLMap(clusters), {}, std::vector<int>(nSymbols + 1, 0), {}};
            for (int i = 0; i < (int) params.hint.size(); ++i) {
                if (params.hint[i] != -1) {
                    if (frand() > 0.5) {
                        hcur.clMap[clusters[i]] = params.hint[i];
                    }
                }
            }
            translate(hcur.clMap, clusters, hcur.plain);
            hcur.memo.resize(N, 1.0);
            hcur.p = calcScore(params, freqMap, hcur.plain, hcur.memo);
            ++nCur;
        }

        // sorted clusters by frequency
        std::vector<std::pair<TClusterId, std::vector<int>>> sorted;

        {
            std::map<TClusterId, std::vector<int>> counts;
            for (int i = 0; i < (int) clusters.size(); ++i) {
                counts[clusters[i]].push_back(i);
            }

            for (auto & kv : counts) {
                sorted.push_back(std::make_pair(kv.first, kv.second));
            }
            std::sort(sorted.begin(), sorted.end(), [](const auto & a, const auto & b) {
                return a.second.size() > b.second.size();
            });
        }

        for (int i = 0; i < (int) sorted.size(); ++i) {
            auto & kvSorted = sorted[i];
            if (kvSorted.second.empty()) break;

            //printf("Processing cluster %2d ('%c') - count = %d\n", kvSorted.first, getEncodedChar(kvSorted.first), (int) kvSorted.second.size());

            if (hypothesesCur[0].clMap.at(kvSorted.first) != 0) {
                continue;
            }

            int nNew = 0;
            for (int j = 0; j < nCur; ++j) {
                const auto & hcur = hypothesesCur[j];

                for (int a = 1; a <= nSymbols; ++a) {
                    // TODO: maybe become parameter
                    // how many clusters can map to the same symbol
                    if (hcur.nused[a] > 20) continue;

                    auto & hnew = hypothesesNew[nNew++];

                    hnew = hcur;
                    hnew.nused[a]++;
                    for (int k = 0; k < (int) kvSorted.second.size(); ++k) {
                        const auto idx = kvSorted.second[k];
                        hnew.plain[idx] = a;

                        const auto idx0 = std::max(0, idx - freqMap.len + 1);
                        const auto idx1 = std::min(N - 1, idx + freqMap.len - 1);
                        std::fill(hnew.memo.begin() + idx0, hnew.memo.begin() + idx1 + 1, 1.0);
                    }
                    hnew.clMap[kvSorted.first] = a;
                    hnew.p = calcScore(params, freqMap, hnew.plain, hnew.memo);
                }
            }

            // sort hypothesesNew by p
            {
                std::vector<std::pair<int, TProb>> sortedNew;
                for (int j = 0; j < nNew; ++j) {
                    sortedNew.push_back(std::make_pair(j, hypothesesNew[j].p));
                }
                std::sort(sortedNew.begin(), sortedNew.end(), [](const auto & a, const auto & b) {
                    return a.second > b.second;
                });

                //printf("Top p = %g\n", (double) hypothesesNew[sortedNew[0].first].p);
                //printDecoded(clusters, hypothesesNew[sortedNew[0].first].clMap);

                nCur = std::min(nHypothesesToKeep, nNew);
                for (int j = 0; j < nCur; ++j) {
                    hypothesesCur[j] = std::move(hypothesesNew[sortedNew[j].first]);
                }
            }
        }

        result.clMap = hypothesesCur[0].clMap;
        result.p = hypothesesCur[0].p;

        return true;
    }

    bool refineNearby(
        const TParameters & params,
        const TFreqMap & freqMap,
        TResult & result) {

        const int N = result.clusters.size();

        struct THypothesis {
            TProb p;
            TClusterToLetterMap clMap;
            std::vector<TLetter> plain;
            std::vector<TProb> memo;
        };

        {
            THypothesis hcur;
            hcur = { 0.0, result.clMap, {}, {}};
            translate(hcur.clMap, result.clusters, hcur.plain);
            for (int i = 0; i < (int) params.hint.size(); ++i) {
                if (params.hint[i] != -1) {
                    hcur.plain[i] = params.hint[i];
                }
            }
            hcur.memo.resize(N, 1.0);
            hcur.p = calcScore(params, freqMap, hcur.plain, hcur.memo);

            THypothesis hbest = hcur;
            while (true) {
                bool converged = true;

                hcur = hbest;
                for (int i = 0; i < N; ++i) {
                    if (hcur.plain[i] < 1 || hcur.plain[i] > 26) continue;
                    auto hnew = hcur;
                    for (int k = 1; k < (int) kNearbyKeys.at('a' + hcur.plain[i] - 1).size(); ++k) {
                        hnew.plain[i] = kNearbyKeys.at('a' + hcur.plain[i] - 1)[k] - 'a' + 1;
                        const auto idx0 = std::max(0, i - freqMap.len + 1);
                        const auto idx1 = std::min(N - 1, i + freqMap.len - 1);
                        std::fill(hnew.memo.begin() + idx0, hnew.memo.begin() + idx1 + 1, 1.0);
                        hnew.p = calcScore(params, freqMap, hnew.plain, hnew.memo);
                        if (hbest.p < hnew.p) {
                            hbest = hnew;
                            converged = false;
                        }
                    }
                }

                if (converged) break;
            }

            printPlain(hbest.plain);
            printf(" [%8.3f %8.3f]", (double) hbest.p, (double) result.p);
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

        //printf("nClusters = %d\n", nClusters);
        for (auto & cid : clusters) {
            assert(cid >= 0 && cid < params.maxClusters);
        }

        return true;
    }

    bool mutateClusters(const TParameters & params, TClusters & clusters) {
        int n = clusters.size();

        for (int i = 0; i < 1; ++i) {
            int idx = rand()%n;

            auto old = clusters[idx];
            do {
                clusters[idx] = 1 + rand()%(params.maxClusters - 1);
            } while (clusters[idx] == old);
        }

        return true;
    }

    double calcPClusters(
            const TParameters & ,
            const TSimilarityMap & ,
            const TSimilarityMap & logMap,
            const TSimilarityMap & logMapInv,
            const TClusters & clusters,
            const TClusterToLetterMap & ) {

        double res = 0.0;
        int n = clusters.size();

        for (int j = 0; j < n - 1; ++j) {
            for (int i = j + 1; i < n; ++i) {
                if (clusters[i] == clusters[j]) {
                    res += logMap[j][i].cc;
                } else {
                    res += logMapInv[j][i].cc;
                }
            }
        }

        return res/((n*(n-1))/2.0);
    }

    bool normalizeSimilarityMap(
            const TParameters & params,
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

        //printf("ccMax = %g, ccMin = %g\n", ccMax, ccMin);

        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                if (i == j) {
                    ccMap[j][i].cc = 1.0;
                    continue;
                }

                auto & v = ccMap[j][i].cc;
                v = (v - ccMin)/(ccMax - ccMin);
                v = std::pow(v, params.fSpread);
                //if (v < 0.50) {
                //    v = 1e-6;
                //}
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

    char getEncodedChar(TClusterId cid) {
        if (cid >= 1 && cid <= 26) {
            return 'a' + cid - 1;
        } else if (cid == 27) {
            return '_';
        }

        return 'A' + cid - 28;
    }

    TLetter decode(const TClusters & t, int idx, const TClusterToLetterMap & clMap, const THint & hint) {
        const auto & cid = t[idx];
        auto let = clMap.at(cid);
        if ((int) hint.size() > idx && hint[idx] != -1) {
            let = hint[idx];
        }

        return let;
    }

    void printEncoded(const TClusters & t) {
        for (const auto & cid : t) {
            printf("%c", getEncodedChar(cid));
        }
    }

    void printDecoded(const TClusters & t, const TClusterToLetterMap & clMap, const THint & hint) {
        for (int i = 0; i < (int) t.size(); ++i) {
            const auto let = decode(t, i, clMap, hint);

            if (let >= 1 && let <= 26) {
                printf("%c", 'a' + let - 1);
            } else if (let == 27) {
                printf("_");
            } else {
                printf(".");
            }
        }
    }

    void printPlain(const std::vector<TLetter> & t) {
        for (int i = 0; i < (int) t.size(); ++i) {
            const auto let = t[i];

            if (let >= 1 && let <= 26) {
                printf("%c", 'a' + let - 1);
            } else if (let == 27) {
                printf("_");
            } else {
                printf(".");
            }
        }
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
        m_curResult = {};

        normalizeSimilarityMap(m_params, m_similarityMap, m_logMap, m_logMapInv);
        generateClustersInitialGuess(m_params, m_similarityMap, m_curResult.clusters);

        //Cipher::beamSearch(m_params, *m_freqMap, m_curResult);
        m_nInitialIters = 0;
        m_pCur = calcPClusters(m_params, m_similarityMap, m_logMap, m_logMapInv, m_curResult.clusters, m_curResult.clMap);
        m_curResult.pClusters = m_pCur;
        m_pZero = m_pCur;

        return true;
    }

    bool Processor::setHint(const THint & hint) {
        m_params.hint = hint;

        return true;
    }

    std::vector<TResult> Processor::getClusterings(int nClusterings) {
        const auto p0 = m_curResult.pClusters;
        //printf("    [getClusterings] p0 = %g\n", p0);

        int nNoImprovement = 0;
        int nTotalIterations = 0;
        auto clustersNew = m_curResult.clusters;
        const int n = clustersNew.size();

        std::vector<TResult> all;
        all.push_back(m_curResult);

        // simulated annealing
        double T = m_params.temp0;
        const double TMin = 0.000001;
        const double alpha = m_params.coolingRate;
        const double pScale = ((n*(n-1))/2.0);

        while (true) {
            clustersNew = m_curResult.clusters;

            // mutate
            int idxChanged = -1;
            {
                idxChanged = rand()%n;

                auto old = clustersNew[idxChanged];
                do {
                    clustersNew[idxChanged] = 1 + rand()%(m_params.maxClusters - 1);
                } while (clustersNew[idxChanged] == old);
            }

            // compute pNew
            auto pNew = m_pCur*pScale;
            {
                const int j = idxChanged;
                for (int i = 0; i < n; ++i) {
                    if (i == j) {
                        continue;
                    }

                    if (m_curResult.clusters[i] == m_curResult.clusters[j]) {
                        pNew -= m_logMap[j][i].cc;
                    } else {
                        pNew -= m_logMapInv[j][i].cc;
                    }

                    if (clustersNew[i] == clustersNew[j]) {
                        pNew += m_logMap[j][i].cc;
                    } else {
                        pNew += m_logMapInv[j][i].cc;
                    }
                }

                pNew /= pScale;
            }

            // check if we should accept the new value
            if (pNew >= m_pCur) {
                m_curResult.clusters = clustersNew;
                m_curResult.pClusters = pNew;
                m_pCur = pNew;
            } else {
                // accept with probability
                const auto pAccept = std::exp((pNew - m_pCur)/T);
                if (pAccept > frand()) {
                    //printf("    [getClusterings] N = %5d, T = %8.8f, pNew = %g, pCur = %g, pAccept = %g\n", nNoImprovement, T, pNew, m_pCur, pAccept);
                    m_curResult.clusters = clustersNew;
                    m_curResult.pClusters = pNew;
                    m_pCur = pNew;
                }
            }

            // check if we should stop
            if (m_pCur > all.back().pClusters) {
                all.push_back(m_curResult);
                nNoImprovement = 0;
            } else {
                nNoImprovement += 1;
            }

            // update temperature
            nTotalIterations++;
            if (nTotalIterations % 1000 == 0) {
                T = T * alpha;
                if (T < TMin) {
                    T = TMin;
                }
                //printf("    [getClusterings] T = %g\n", T);
            }

            if (nNoImprovement > 1000 && T < 2*TMin) {
                break;
            }
        }

        //printf("    [getClusterings] nTotalIterations = %d\n", nTotalIterations);
        //printf("    [getClusterings] pFinal = %g\n", m_curResult.pClusters);

        const auto pClustersBest = all.back().pClusters;

        std::vector<TResult> result;
        {
            result.push_back(all.back());

            for (int i = 1; i < nClusterings; ++i) {
                int jBest = -1;
                double pDiffMax = -1.0;
                for (int j = 5*all.size()/6; j < (int) all.size(); ++j) {
                    if (all[j].pClusters < 1.1*pClustersBest) {
                        continue;
                    }
                    double pDiffMin = std::numeric_limits<double>::max();
                    for (int k = 0; k < (int) result.size(); ++k) {
                        double pDiff = std::fabs(all[j].pClusters - result[k].pClusters);
                        if (pDiffMin > pDiff) {
                            pDiffMin = pDiff;
                        }
                    }
                    if (pDiffMax < pDiffMin) {
                        pDiffMax = pDiffMin;
                        jBest = j;
                    }
                }

                if (pDiffMax < 0.005 || jBest == -1) {
                    break;
                }

                result.push_back(all[jBest]);
            }
        }

        std::sort(result.begin(), result.end(), [](const TResult & a, const TResult & b) {
            return a.pClusters > b.pClusters;
        });

        return result;
    }

    bool Processor::compute() {
        auto clusterings = getClusterings(1);
        m_curResult = clusterings[0];
        Cipher::beamSearch(m_params, *m_freqMap, m_curResult);
        m_curResult.id++;

        return true;
    }

    const TResult & Processor::getResult() const {
        return m_curResult;
    }

    const TSimilarityMap & Processor::getSimilarityMap() const {
        return m_similarityMap;
    }

    float findBestCutoffFreq(const TWaveformF & waveform, EAudioFilter filterId, int64_t sampleRate, float minCutoffFreq_Hz, float maxCutoffFreq_Hz, float step_Hz) {
        double pClustersBest = -1e10;
        float freqCutoffBest_Hz = minCutoffFreq_Hz;

        for (float freqCutoff_Hz = minCutoffFreq_Hz; freqCutoff_Hz <= maxCutoffFreq_Hz; freqCutoff_Hz += step_Hz) {
            TWaveformI16 waveformInput;
            TWaveformF waveformFiltered = waveform;
            filter(waveformFiltered, (EAudioFilter) filterId, freqCutoff_Hz, sampleRate);

            if (convert(waveformFiltered, waveformInput) == false) {
                fprintf(stderr, "%s:%d: convert() failed\n", __FILE__, __LINE__);
                return minCutoffFreq_Hz;
            }

            TKeyPressCollectionI16 keyPresses;
            {
                TWaveformI16 waveformMax;
                TWaveformI16 waveformThreshold;
                if (findKeyPresses(getView(waveformInput, 0), keyPresses, waveformThreshold, waveformMax,
                                   kFindKeysThreshold, kFindKeysHistorySize, kFindKeysHistorySizeReset, kFindKeysRemoveLowPower) == false) {
                    fprintf(stderr, "%s:%d: findKeyPresses() failed\n", __FILE__, __LINE__);
                    return minCutoffFreq_Hz;
                }
            }

            TSimilarityMap similarityMap;
            if (calculateSimilartyMap(kKeyWidth_samples, kKeyAlign_samples, kKeyWidth_samples - kKeyOffset_samples, keyPresses, similarityMap) == false) {
                fprintf(stderr, "%s:%d: calculateSimilartyMap() failed\n", __FILE__, __LINE__);
                return minCutoffFreq_Hz;
            }

            {
                Cipher::TFreqMap freqMap; // not used for anything
                Cipher::Processor processor;

                Cipher::TParameters params;
                params.maxClusters = 50;
                params.wEnglishFreq = 20.0;
                processor.init(params, freqMap, similarityMap);

                auto clusteringsCur = processor.getClusterings(1);
                if (clusteringsCur[0].pClusters > pClustersBest) {
                    pClustersBest = clusteringsCur[0].pClusters;
                    freqCutoffBest_Hz = freqCutoff_Hz;
                }
                printf("    [findBestCutoffFreq] freqCutoff_Hz = %g, pClusters = %g\n", freqCutoff_Hz, clusteringsCur[0].pClusters);
            }
        }

        return freqCutoffBest_Hz;
    }

}
