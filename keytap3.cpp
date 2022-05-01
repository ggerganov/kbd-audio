/*! \file keytap3.cpp
 *  \brief Fully automated acoustic keyboard eavesdropping
 *  \author Georgi Gerganov
 */

#include "common.h"
#include "constants.h"
#include "subbreak3.h"

#include <chrono>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include <mutex>
#include <thread>

#define MY_DEBUG

using TSampleInput          = TSampleF;
using TSample               = TSampleI16;
using TWaveform             = TWaveformI16;
using TWaveformView         = TWaveformViewI16;
using TKeyPressData         = TKeyPressDataI16;
using TKeyPressCollection   = TKeyPressCollectionI16;

int main(int argc, char ** argv) {
    printf("Usage: %s record.kbd n-gram-dir [-FN] [-fN]\n", argv[0]);
    printf("    -FN - select filter type, (0 - none, 1 - first order high-pass, 2 - second order high-pass)\n");
    printf("    -fN - cutoff frequency in Hz\n");
    if (argc < 3) {
        return -1;
    }

    const int64_t sampleRate = kSampleRate;

    const auto argm = parseCmdArguments(argc, argv);
    const int filterId      = argm.count("F") == 0 ? EAudioFilter::FirstOrderHighPass : std::stoi(argm.at("F"));

    int freqCutoff_Hz = argm.count("f") == 0 ? 0 : std::stoi(argm.at("f"));

    Cipher::TFreqMap freqMap6;
    {
        const auto tStart = std::chrono::high_resolution_clock::now();

        printf("[+] Loading n-grams from '%s'\n", argv[2]);

        if (Cipher::loadFreqMapBinary((std::string(argv[2]) + "/ggwords-6-gram.dat.binary").c_str(), freqMap6) == false) {
            return -5;
        }

        const auto tEnd = std::chrono::high_resolution_clock::now();

        printf("[+] Loading took %4.3f seconds\n", toSeconds(tStart, tEnd));
    }

    // Main algorithm

    TWaveform waveformInput;
    {
        TWaveformF waveformInputF;
        printf("[+] Loading recording from '%s'\n", argv[1]);
        if (readFromFile<TSampleF>(argv[1], waveformInputF) == false) {
            printf("Specified file '%s' does not exist\n", argv[1]);
            return -1;
        } else {
            if (freqCutoff_Hz == 0) {
                const auto tStart = std::chrono::high_resolution_clock::now();

                freqCutoff_Hz = Cipher::findBestCutoffFreq(waveformInputF, (EAudioFilter) filterId, kSampleRate, 100.0f, 1000.0f, 100.0f);

                const auto tEnd = std::chrono::high_resolution_clock::now();
                printf("[+] Found best freqCutoff = %d Hz, took %4.3f seconds\n", freqCutoff_Hz, toSeconds(tStart, tEnd));
            }

            printf("[+] Filtering waveform with filter type = %d and cutoff frequency = %d Hz\n", filterId, freqCutoff_Hz);
            ::filter(waveformInputF, (EAudioFilter) filterId, freqCutoff_Hz, kSampleRate);

            printf("[+] Converting waveform to i16 format ...\n");
            if (convert(waveformInputF, waveformInput) == false) {
                printf("Conversion failed\n");
                return -4;
            }
        }
    }

    printf("[+] Loaded recording: of %d samples (sample size = %d bytes)\n", (int) waveformInput.size(), (int) sizeof(TSample));
    printf("    Size in memory:          %g MB\n", (float)(sizeof(TSample)*waveformInput.size())/1024/1024);
    printf("    Sample size:             %d\n", (int) sizeof(TSample));
    printf("    Total number of samples: %d\n", (int) waveformInput.size());
    printf("    Recording length:        %g seconds\n", (float)(waveformInput.size())/sampleRate);

    TKeyPressCollection keyPresses;
    {
        const auto tStart = std::chrono::high_resolution_clock::now();

        printf("[+] Searching for key presses\n");

        TWaveform waveformMax;
        TWaveform waveformThreshold;
        if (findKeyPresses(getView(waveformInput, 0), keyPresses, waveformThreshold, waveformMax,
                           kFindKeysThreshold, kFindKeysHistorySize, kFindKeysHistorySizeReset, kFindKeysRemoveLowPower) == false) {
            printf("Failed to detect keypresses\n");
            return -2;
        }

        const auto tEnd = std::chrono::high_resolution_clock::now();

        printf("[+] Detected a total of %d potential key presses\n", (int) keyPresses.size());
        printf("[+] Search took %4.3f seconds\n", toSeconds(tStart, tEnd));
    }

    int n = keyPresses.size();

    TSimilarityMap similarityMap;
    {
        const auto tStart = std::chrono::high_resolution_clock::now();

        printf("[+] Calculating CC similarity map\n");

        if (calculateSimilartyMap(kKeyWidth_samples, kKeyAlign_samples, kKeyWidth_samples - kKeyOffset_samples, keyPresses, similarityMap) == false) {
            printf("Failed to calculate similariy map\n");
            return -3;
        }

        const auto tEnd = std::chrono::high_resolution_clock::now();

        printf("[+] Calculation took %4.3f seconds\n", toSeconds(tStart, tEnd));

        {
            const auto tStart = std::chrono::high_resolution_clock::now();

            printf("[+] Removing low-similarity keys\n");

            const int n0 = keyPresses.size();

            if (removeLowSimilarityKeys(keyPresses, similarityMap, 0.3f) == false) {
                printf("Failed to remove low-similarity keys\n");
                return -4;
            }

            const int n1 = keyPresses.size();

            const auto tEnd = std::chrono::high_resolution_clock::now();

            printf("[+] Removed %d low-similarity keys, took %4.3f seconds\n", n0 - n1, toSeconds(tStart, tEnd));
        }

        n = keyPresses.size();

        const int ncc = std::min(32, n);
        for (int j = 0; j < ncc; ++j) {
            printf("%2d: ", j);
            for (int i = 0; i < ncc; ++i) {
                printf("%6.3f ", similarityMap[j][i].cc);
            }
            printf("\n");
        }
        printf("\n");

        auto minCC = similarityMap[0][1].cc;
        auto maxCC = similarityMap[0][1].cc;
        for (int j = 0; j < n - 1; ++j) {
            for (int i = j + 1; i < n; ++i) {
                minCC = std::min(minCC, similarityMap[j][i].cc);
                maxCC = std::max(maxCC, similarityMap[j][i].cc);
            }
        }

        printf("[+] Similarity map: min = %g, max = %g\n", minCC, maxCC);
    }

    printf("[+] Attempting to recover the text from the recording ...\n");

    for (int iMain = 0; iMain < 16; ++iMain) {
        Cipher::Processor processor;

        Cipher::TParameters params;
        params.maxClusters = 30;
        params.wEnglishFreq = 30.0;
        params.fSpread = 0.5 + 0.1*iMain;
        params.nHypothesesToKeep = std::max(100, 500 - 2*std::min(200, std::max(0, ((int) keyPresses.size() - 100))));
        processor.init(params, freqMap6, similarityMap);

        std::vector<Cipher::TResult> clusterings;

        // clustering
        {
            const auto tStart = std::chrono::high_resolution_clock::now();

            for (int nIter = 0; nIter < 16; ++nIter) {
                auto clusteringsCur = processor.getClusterings(2);

                for (int i = 0; i < (int) clusteringsCur.size(); ++i) {
                    clusterings.push_back(std::move(clusteringsCur[i]));
                }

                params.maxClusters = 30 + 4*(nIter + 1);
                processor.init(params, freqMap6, similarityMap);
            }

            const auto tEnd = std::chrono::high_resolution_clock::now();
            printf("[+] Clustering took %4.3f seconds\n", toSeconds(tStart, tEnd));
        }

        params.hint.clear();
        params.hint.resize(n, -1);

        // beam search
        int nThread = std::min((int) std::thread::hardware_concurrency(), (int) clusterings.size());
        {
            std::vector<std::thread> workers(nThread);

            std::mutex mutexPrint;
            for (int i = 0; i < nThread; ++i) {
                workers[i] = std::thread([&, i]() {
                    for (int j = i; j < (int) clusterings.size(); j += nThread) {
                        Cipher::beamSearch(params, freqMap6, clusterings[j]);
                        mutexPrint.lock();
                        printf(" ");
                        Cipher::printDecoded(clusterings[j].clusters, clusterings[j].clMap, params.hint);
                        printf(" [%8.3f %8.3f]\n", clusterings[j].p, clusterings[j].pClusters);
                        mutexPrint.unlock();
                    }
                });
            }

            for (auto& worker : workers) {
                worker.join();
            }
        }
    }

    return 0;
}
