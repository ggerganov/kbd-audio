/*! \file subbreak2.h
 *  \brief Breaking non-exact substitution ciphers
 *  \author Georgi Gerganov
 */

#pragma once

#include "common.h"

#include <map>
#include <cmath>
#include <vector>
#include <string>

namespace Cipher {

    using TCode = int32_t;
    using TProb = double;
    using TGramLen = int32_t;
    using THint = std::vector<int32_t>;
    using TPlainText = std::vector<TLetter>;
    using TLetterCount = std::vector<int>;
    using TClusterCount = std::vector<int>;
    using TClusterPos = std::vector<std::vector<int>>;

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
        int nSubbreakIterations = 300;
        TProb pNonAlphabetic = 0.0;

        // Metropolis–Hastings params
        int nMHInitialIters = 100;
        int nMHIters = 100;

        // language model
        bool includeSpaces = true;
        float wEnglishFreq = 10.0f;
        float wLanguageModel = 1.0f;

        // simulation parameters
        float waveformDeviationMin = 0.1f;
        float waveformDeviationSig = 0.2f;
        float waveformDetectionErrorP = 0.10f;
        float waveformDetectionErrorMin = 0.5f;
        float waveformDetectionErrorSig = 0.4f;
        float similarityNoiseSig = 0.1f;
        float similarityMismatchAvg = 0.3f;
        float similarityMismatchSig = 0.2f;

        THint hint = {};
    };

    struct TFreqMap {
        TGramLen len = -1;
        TProb pmin = 0.0;
        int64_t nTotal = 0;
        std::vector<TProb> prob;
    };

    struct TResult {
        int32_t id = 0;
        TProb p = -100.0;
        TClusterToLetterMap clMap;
        TClusters clusters;
    };

    TCode calcCode(const char * data, int n);
    bool loadFreqMap(const char * fname, TFreqMap & res);

    bool encryptExact(const TParameters & params, const std::string & text, TClusters & clusters);
    bool generateSimilarityMap(const TParameters & params, const std::string & text, TSimilarityMap & ccMap);

    bool generateClusters(const TParameters & params, int n, TClusters & clusters);
    bool printClusterGoodness(const std::string & text, const TClusters & clusters);

    float costF(const TSimilarityMap & ccMap, const TClusters & clusters);
    float costFUpdate(const TSimilarityMap & ccMap, const TClusters & clusters, int i, int cid, float c0);

    bool doSimulatedAnnealing3(
            const TParameters & params,
            const TSimilarityMap & ccMap,
            TClusters & clusters);

    bool doSimulatedAnnealing4(
            const TParameters & params,
            const TFreqMap & freqMap,
            const TClusters & clusters,
            TClusterToLetterMap & clMap);

    bool doSimulatedAnnealing5(
            const TParameters & params,
            const TFreqMap & freqMap,
            const TSimilarityMap & ccMap,
            TClusters & clusters,
            TClusterToLetterMap & clMap);

    bool subbreak(
            const TParameters & params,
            const TFreqMap & freqMap,
            TResult & result);

    bool subbreak1(
            const TParameters & params,
            const TFreqMap & freqMap,
            TResult & result);

    bool generateClustersInitialGuess(
            const TParameters & params,
            const TSimilarityMap & ccMap,
            TClusters & clusters);

    bool mutateClusters(const TParameters & params, TClusters & clusters);

    double calcPClusters(
            const TParameters & ,
            const TSimilarityMap & ,
            const TSimilarityMap & logMap,
            const TSimilarityMap & logMapInv,
            const TClusters & clusters,
            const TClusterToLetterMap & clMap);

    bool normalizeSimilarityMap(
            const TParameters & ,
            TSimilarityMap & ccMap,
            TSimilarityMap & logMap,
            TSimilarityMap & logMapInv);

    void printText(const TClusters & t);
    void printText(const TClusters & t, const TClusterToLetterMap & clMap);

    class Processor {
    public:
        Processor();

        bool init(
                const TParameters & params,
                const TFreqMap & freqMap,
                const TSimilarityMap & similarityMap);

        bool setPNonAlphabetic(TProb p);
        bool setWEnglishFreq(float w);
        bool setHint(const THint & hint);

        bool compute();

        int getIters() const { return m_nMHInitialIters; }
        const TResult & getResult() const;
        const TSimilarityMap & getSimilarityMap() const;

    private:
        TParameters m_params;
        const TFreqMap* m_freqMap = nullptr;
        TSimilarityMap m_similarityMap;
        TSimilarityMap m_logMap;
        TSimilarityMap m_logMapInv;

        int m_nMHInitialIters = 0;
        double m_pCur = 0.0f;

        TResult m_curResult;
    };
}
