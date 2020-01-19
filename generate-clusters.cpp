/*! \file generate-clusters.cpp
 *  \brief Generate clusters given a similarity matrix using metropolis-hastings algorithm
 *  \author Georgi Gerganov
 */

#include "subbreak2.h"

#include <cassert>
#include <random>
#include <algorithm>

inline float frand() { return ((float)rand())/RAND_MAX; }

inline float frandGaussian(float mu, float sigma)
{
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

void generateClustersRandom(const Cipher::TParameters & params, int n, Cipher::TClusters & clusters) {
    clusters.resize(n);
    for (int i = 0; i < n; ++i) {
        clusters[i] = rand()%params.maxClusters;
    }
}

void generateClustersInitialGuess(
        const Cipher::TParameters & params,
        const Cipher::TSimilarityMap & ccMap,
        Cipher::TClusters & clusters) {
    int n = ccMap.size();

    int nClusters = n;
    clusters.resize(n);
    for (int i = 0; i < n; ++i) {
        clusters[i] = i;
    }

    struct Pair {
        int i;
        int j;
        float cc;

        bool operator < (const Pair & a) const { return cc > a.cc; }
    };

    std::vector<Pair> ccPairs;
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            ccPairs.emplace_back(Pair{i, j, ccMap.at(i).at(j)});
        }
    }

    std::sort(ccPairs.begin(), ccPairs.end());

    for (int k = 0; k < (int) ccPairs.size(); ++k) {
        int i = ccPairs[k].i;
        int j = ccPairs[k].j;

        if (clusters[i] == clusters[j]) continue;

        int cidi = clusters[i];
        int cidj = clusters[j];
        if (cidi > cidj) {
            std::swap(cidi, cidj);
        }

        for (int p = 0; p < n; ++p) {
            if (clusters[p] == cidj) {
                clusters[p] = cidi;
            }
        }
        --nClusters;

        if (nClusters <= params.maxClusters) break;
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
}

void mutateClusters(const Cipher::TParameters & params, Cipher::TClusters & clusters) {
    int n = clusters.size();
    int k = std::max(1, std::min(n, (int) std::fabs(std::round(frandGaussian(0.0, 3)))));

    auto midxs = subset(k, n);

    for (auto & idx : midxs) {
        auto p = clusters[idx];
        while (clusters[idx] == p) {
            clusters[idx] = rand()%params.maxClusters;
        }
    }
}

float calcPClusters(
        const Cipher::TParameters & params,
        const Cipher::TSimilarityMap & ccMap,
        const Cipher::TSimilarityMap & logMap,
        const Cipher::TSimilarityMap & logMapInv,
        const Cipher::TClusters & clusters) {
    float res = 0.0;
    int n = clusters.size();

    for (int j = 0; j < n - 1; ++j) {
        for (int i = j + 1; i < n; ++i) {
            if (clusters[i] == clusters[j]) {
                res += logMap.at(j).at(i);
            } else {
                res += logMapInv.at(j).at(i);
            }

            //if (clusters[i] == clusters[j]) {
            //    res += ccMap.at(j).at(i);
            //} else {
            //    res += 1.0f - ccMap.at(j).at(i);
            //}
        }
    }

    //res /= (n*(n-1))/2;

    return res;
}

void normalizeSimilarityMap(
        const Cipher::TParameters & params,
        Cipher::TSimilarityMap & ccMap,
        Cipher::TSimilarityMap & logMap,
        Cipher::TSimilarityMap & logMapInv) {
    int n = ccMap.size();

    float ccMin = std::numeric_limits<float>::max();
    float ccMax = std::numeric_limits<float>::min();

    for (int j = 0; j < n - 1; ++j) {
        for (int i = j + 1; i < n; ++i) {
            ccMin = std::min(ccMin, ccMap[j][i]);
            ccMax = std::max(ccMax, ccMap[j][i]);
        }
    }

    ccMin -= 1e-6;
    ccMax += 1e-6;

    printf("ccMax = %g, ccMin = %g\n", ccMax, ccMin);

    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            if (i == j) {
                ccMap[j][i] = 1.0f;
                continue;
            }

            auto & v = ccMap[j][i];
            v = (v - ccMin)/(ccMax - ccMin);
            //v = 1.0 - std::exp(-1.1f*v);
        }
    }

    logMap = ccMap;

    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            if (i == j) {
                logMap[j][i] = 0.0f;
                continue;
            }

            auto & v = logMap[j][i];
            v = std::log(v);
        }
    }

    logMapInv = ccMap;

    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            if (i == j) {
                logMap[j][i] = -1e6;
                continue;
            }

            auto & v = logMapInv[j][i];
            v = std::log(1.0f - v);
        }
    }
}

int main(int argc, char ** argv) {
    printf("Usage: %s n-gram.txt\n", argv[0]);
    if (argc < 2) {
        return -1;
    }

    //srand(time(0));

    Cipher::TFreqMap freqMap;
    if (Cipher::loadFreqMap(argv[1], freqMap) == false) {
        return -1;
    }

    std::string plain;

    plain = R"(
as far as services go, only two steady contributors of revenue
streams keep swelling without apple having to charge
subscription fees. one is the money paid by google parent alphabet
for searches made through apple products such as the safari browser and siri.
    )";

    Cipher::TParameters params;
    Cipher::TSimilarityMap ccMap;

    Cipher::generateSimilarityMap(params, plain, ccMap);

    Cipher::TSimilarityMap logMap;
    Cipher::TSimilarityMap logMapInv;
    normalizeSimilarityMap(params, ccMap, logMap, logMapInv);

    Cipher::TClusters clusters;
    //generateClustersRandom(params, plain.size(), clusters);
    generateClustersInitialGuess(params, ccMap, clusters);

    auto pCur = calcPClusters(params, ccMap, logMap, logMapInv, clusters);
    while (true) {
        auto clustersNew = clusters;
        mutateClusters(params, clustersNew);
        auto pNew = calcPClusters(params, ccMap, logMap, logMapInv, clustersNew);

        //printf("pNew = %g, pCur = %g\n", pNew, pCur);

        auto u = frand();
        //auto alpha = pNew/pCur;
        auto alpha = std::exp((pNew - pCur));

        //printf("alpha = %g\n", alpha);

        if (u <= alpha) {
            clusters = clustersNew;
            pCur = pNew;

            int n = plain.size();
            int nMatch = 0;
            for (int j = 0; j < n - 1; ++j) {
                for (int i = j + 1; i < n; ++i) {
                    if ((plain[i] == plain[j] && clusters[i] == clusters[j]) ||
                        (plain[i] != plain[j] && clusters[i] != clusters[j])) {
                        ++nMatch;
                    }
                }
            }

            float pMatch = float(nMatch)/((n*(n-1))/2);

            printf("pCur = %g, pMatch = %g, alpha = %g\n", pCur, pMatch, alpha);

            static int cnt = 99;
            if (++cnt >= 100 && pMatch > 0.94) {
                cnt = 0;
                params.nSubbreakIterations = 1000;
                Cipher::TClusterToLetterMap clMap;
                Cipher::subbreak(params, freqMap, clusters, clMap);
            }
        }
    }

    return 0;
}
