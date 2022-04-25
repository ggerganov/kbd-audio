/*! \file keytap3-multi.cpp
 *  \brief Same as keytap3.cpp but using multi-samples
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

#define MY_DEBUG

using TSampleInput          = TSampleF;
using TSample               = TSampleI16;

int main(int argc, char ** argv) {
    srand(time(0));

    printf("Usage: %s record.kbd n-gram-dir [-FN] [-fN]\n", argv[0]);
    printf("    -FN - select filter type, (0 - none, 1 - first order high-pass, 2 - second order high-pass)\n");
    printf("    -fN - cutoff frequency in Hz\n");
    if (argc < 3) {
        return -1;
    }

    const int64_t sampleRate = kSampleRate;

    const auto argm = parseCmdArguments(argc, argv);
    const int filterId      = argm.count("F") == EAudioFilter::FirstOrderHighPass ? 0 : std::stoi(argm.at("F"));
    const int freqCutoff_Hz = argm.count("f") == kFreqCutoff_Hz ? 0 : std::stoi(argm.at("f"));

    TWaveformMI16 waveformInputMI16;
    {
        TWaveformF waveformInputF;
        printf("[+] Loading recording from '%s'\n", argv[1]);
        if (readFromFile<TSampleF>(argv[1], waveformInputF) == false) {
            printf("Specified file '%s' does not exist\n", argv[1]);
            return -1;
        } else {
            for (int j = 0; j < TSampleMI16::N; ++j) {
                TWaveformI16 waveformInputI16;

                printf("[+] Filtering waveform with filter type = %d and cutoff frequency = %d Hz\n", filterId, freqCutoff_Hz + j*200);
                auto filtered = waveformInputF;
                ::filter(filtered, (EAudioFilter) (j%2 + filterId), freqCutoff_Hz + j*200, kSampleRate);

                //{
                //    auto tmp = filtered;
                //    for (int i = j; i < (int) filtered.size() - j; ++i) {
                //        float sum = 0.0f;
                //        for (int w = -j; w <= j; ++w) {
                //            sum += (tmp[i + w]);
                //        }
                //        filtered[i] = sum / (2 * j + 1);
                //    }
                //}

                printf("[+] Converting waveform to i16 format ...\n");
                if (convert(filtered, waveformInputI16) == false) {
                    printf("Conversion failed\n");
                    return -4;
                }

                waveformInputMI16.resize(waveformInputI16.size());
                for (size_t i = 0; i < waveformInputI16.size(); ++i) {
                    waveformInputMI16[i][j] = waveformInputI16[i];
                }
            }
        }
    }

    printf("[+] Loaded recording: of %d samples (sample size = %d bytes)\n", (int) waveformInputMI16.size(), (int) sizeof(TSample));
    printf("    Size in memory:          %g MB\n", (float)(sizeof(TSample)*waveformInputMI16.size())/1024/1024);
    printf("    Sample size:             %d\n", (int) sizeof(TSample));
    printf("    Total number of samples: %d\n", (int) waveformInputMI16.size());
    printf("    Recording length:        %g seconds\n", (float)(waveformInputMI16.size())/sampleRate);

    TKeyPressCollectionMI16 keyPresses;
    {
        const auto tStart = std::chrono::high_resolution_clock::now();

        printf("[+] Searching for key presses\n");

        TWaveformMI16 waveformMax;
        TWaveformMI16 waveformThreshold;
        if (findKeyPresses(getView(waveformInputMI16, 0), keyPresses, waveformThreshold, waveformMax, 8.0, 512, 2*1024, true) == false) {
            printf("Failed to detect keypresses\n");
            return -2;
        }

        const auto tEnd = std::chrono::high_resolution_clock::now();

        printf("[+] Detected a total of %d potential key presses\n", (int) keyPresses.size());
        printf("[+] Search took %4.3f seconds\n", toSeconds(tStart, tEnd));
    }

    TSimilarityMap similarityMap;
    {
        const auto tStart = std::chrono::high_resolution_clock::now();

        printf("[+] Calculating CC similarity map\n");

        if (calculateSimilartyMap(256, 32, 128, keyPresses, similarityMap) == false) {
            printf("Failed to calculate similariy map\n");
            return -3;
        }

        const auto tEnd = std::chrono::high_resolution_clock::now();

        printf("[+] Calculation took %4.3f seconds\n", toSeconds(tStart, tEnd));

        {
            const auto nPrint = std::min(16, (int) keyPresses.size());
            for (int j = 0; j < nPrint; ++j) {
                printf("%2d: ", j);
                for (int i = 0; i < nPrint; ++i) {
                    printf("%6.3f ", similarityMap[j][i].cc);
                }
                printf("\n");
            }
            printf("\n");
        }
    }

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

    {
        Cipher::Processor processor;

        Cipher::TParameters params;
        params.maxClusters = 40;
        params.nIters = 1;
        params.nInitialIters = 0;
        params.wEnglishFreq = 10.0f;
        processor.init(params, freqMap6, similarityMap);

        double lastP = -1000.0;
        int nNoImprovement = 0;

        printf("[+] Attempting to recover the text from the recording ...\n");

        while (true) {
            processor.compute();

            {
                const auto & result = processor.getResult();
                if (result.p < -100.0f) continue;

                if (lastP < result.p) {
                    lastP = result.p;
                    nNoImprovement = 0;
                    printf("%8.3f %8.3f ", result.p, result.pClusters);
                    Cipher::printDecoded(result.clusters, result.clMap, {});
                } else {
                    nNoImprovement++;
                }

            }

            if (nNoImprovement < 1000) {
                continue;
            }

            //printf("[+] Restarting\n");

            lastP = -1000.0;
            nNoImprovement = 0;

            params.maxClusters = 40 + rand()%50;
            params.nInitialIters = 100 + rand()%5000;
            params.wEnglishFreq = 2.0f + rand()%50;
            processor.init(params, freqMap6, similarityMap);
        }
    }

    return 0;
}
