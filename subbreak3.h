/*! \file subbreak3.h
 *  \brief Breaking non-exact substitution ciphers
 *  \author Georgi Gerganov
 */

#pragma once

#include "common.h"

#include <map>
#include <cmath>
#include <vector>
#include <string>
#include <unordered_map>

namespace Cipher {

    using TCode = int32_t;
    using TProb = float;
    using TGramLen = int32_t;
    using THint = std::vector<int32_t>;
    using TPlainText = std::vector<TLetter>;
    using TLetterCount = std::vector<int>;
    using TClusterCount = std::vector<int>;
    using TClusterPos = std::vector<std::vector<int>>;

    struct TParameters {
        // clustering params
        int minClusters = 23;
        int maxClusters = 32;
        int nInitialIters = 1000;
        int nIters = 100;
        double fSpread = 1.0;

        // simulated annealing params
        float temp0 = 0.0001;
        float coolingRate = 0.95;

        // language model
        float wEnglishFreq = 10.0f;

        // beam search
        int nHypothesesToKeep = 500;

        THint hint = {};
    };

    struct TFreqMap {
        TGramLen len = -1;
        int64_t nTotal = 0;
        TProb pmin = 0;
        std::unordered_map<TCode, TProb> prob;
    };

    struct TResult {
        int32_t id = 0;
        TProb p = -999.0;
        double pClusters = -999.0;
        TClusterToLetterMap clMap;
        TClusters clusters;
    };

    TCode calcCode(const char * data, int n);

    // n-grams with lower probability than pmin are assigned cost = log10(pmin)
    bool loadFreqMap(const char * fname, TFreqMap & res, double pmin = 0.000001);

    bool saveFreqMapBinary(const char * fname, const TFreqMap & res);
    bool loadFreqMapBinary(const char * fname, TFreqMap & res);

    bool encryptExact(const TParameters & params, const std::string & text, TClusters & clusters);

    bool beamSearch(
            const TParameters & params,
            const TFreqMap & freqMap,
            TResult & result);

    bool refineNearby(
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

    char getEncodedChar(TClusterId);

    TLetter decode(const TClusters & t, int idx, const TClusterToLetterMap & clMap, const THint & hint);

    void printEncoded(const TClusters & t);
    void printDecoded(const TClusters & t, const TClusterToLetterMap & clMap, const THint & hint);
    void printPlain(const std::vector<TLetter> & t);

    class Processor {
    public:
        Processor();

        bool init(
                const TParameters & params,
                const TFreqMap & freqMap,
                const TSimilarityMap & similarityMap);

        bool setHint(const THint & hint);

        std::vector<TResult> getClusterings(int nClusterings);
        bool compute();

        int getIters() const { return m_nInitialIters; }
        const TResult & getResult() const;
        const TSimilarityMap & getSimilarityMap() const;

    private:
        TParameters m_params;
        const TFreqMap* m_freqMap = nullptr;
        TSimilarityMap m_similarityMap;
        TSimilarityMap m_logMap;
        TSimilarityMap m_logMapInv;

        int m_nInitialIters = 0;
        double m_pCur = 0.0f;
        double m_pZero = 0.0f;

        TResult m_curResult;
    };

    float findBestCutoffFreq(
            const TWaveformF & waveform,
            EAudioFilter filterId,
            int64_t sampleRate,
            float minCutoffFreq_Hz,
            float maxCutoffFreq_Hz,
            float step_Hz);
}
