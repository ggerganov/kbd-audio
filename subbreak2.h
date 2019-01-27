/*! \file subbreak2.h
 *  \brief Breaking non-exact substituition ciphers
 *  \author Georgi Gerganov
 */

#pragma once

#include <map>
#include <cmath>
#include <vector>
#include <string>

namespace Cipher {

    using TCode = int32_t;
    using TProb = double;
    using TGramLen = int32_t;
    using TLetter = int32_t;
    using TCluster = int32_t;
    using TClusters = std::vector<TCluster>;
    using TSimilarityMap = std::map<int, std::map<int, float>>;
    using TClusterToLetterMap = std::map<TCluster, TLetter>;

    struct TParameters {
        // clustering params
        int minClusters = 23;
        int maxClusters = 30;
        int nIterationsPerUpdate = 100000;
        int nChangePerIteration = 1;

        // simulated annealing params
        int saMaxIterations = 20000000;
        float temp0 = 1e0;
        float coolingRate = 0.999995;

        // subbreak params
        int nSubbreakIterations = 100000;

        // language model
        float wEnglishFreq = 10.0f;
        float wLanguageModel = 1.0f;
    };

    struct TFreqMap {
        TGramLen len;
        TProb pmin;
        std::vector<TProb> prob;
    };

    TCode calcCode(const char * data, int n);
    bool loadFreqMap(const char * fname, TFreqMap & res);

    bool generateClusters(const TParameters & params, int N, TClusters & clusters, const std::vector<int> & hint);
    bool printClusterGoodness(const std::string & text, const TClusters & clusters);

    float costF(const TSimilarityMap & ccMap, const TClusters & clusters);
    float costFUpdate(const TSimilarityMap & ccMap, const TClusters & clusters, int i, int cid, float c0);
    TProb calcScore0(const TFreqMap & freqMap, const TClusters & txt, const TClusterToLetterMap & clMap);

    bool doSimulatedAnnealing3(
        const TParameters & params,
        const TSimilarityMap & ccMap,
        TClusters & clusters,
        const std::vector<int> & hint = {}
        );

    bool doSimulatedAnnealing4(
        const TParameters & params,
        const TFreqMap & freqMap,
        const TClusters & clusters,
        TClusterToLetterMap & clMap,
        const std::vector<int> & hint = {}
        );

    bool subbreak(
        const TParameters & params,
        const TFreqMap & freqMap,
        const TClusters & clusters,
        TClusterToLetterMap & clMap,
        const std::vector<int> & hint = {}
        );

    bool doSimulatedAnnealing5(
        const TParameters & params,
        const TFreqMap & freqMap,
        const TSimilarityMap & ccMap,
        TClusters & clusters,
        TClusterToLetterMap & clMap,
        const std::vector<int> & hint = {}
        );
}
