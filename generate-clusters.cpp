/*! \file generate-clusters.cpp
 *  \brief Generate clusters given a similarity matrix using metropolis-hastings algorithm
 *  \author Georgi Gerganov
 */

#include "subbreak2.h"

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
    TSimilarityMap ccMap;

    Cipher::generateSimilarityMap(params, plain, ccMap);

    TSimilarityMap logMap;
    TSimilarityMap logMapInv;
    Cipher::normalizeSimilarityMap(params, ccMap, logMap, logMapInv);

    TClusters clusters;
    Cipher::generateClustersInitialGuess(params, ccMap, clusters);

    auto pCur = Cipher::calcPClusters(params, ccMap, logMap, logMapInv, clusters);
    while (true) {
        auto clustersNew = clusters;
        Cipher::mutateClusters(params, clustersNew);
        auto pNew = Cipher::calcPClusters(params, ccMap, logMap, logMapInv, clustersNew);

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
                TClusterToLetterMap clMap;
                Cipher::subbreak(params, freqMap, clusters, clMap);
            }
        }
    }

    return 0;
}
