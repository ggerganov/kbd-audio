/*! \file key-detector.cpp
 *  \brief Capture audio until N keystrokes are detected
 *  \author Georgi Gerganov
 */

#include "constants.h"
#include "common.h"
#include "audio-logger.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <chrono>
#include <thread>
#include <fstream>

int main(int argc, char ** argv) {
    printf("Usage: %s output.kbd nkeys [-cN] [-CN]\n", argv[0]);
    printf("    -cN - select capture device N\n");
    printf("    -CN - number N of capture channels N\n");
    printf("\n");

    if (argc < 3) {
        return -127;
    }

    const auto argm = parseCmdArguments(argc, argv);
    const int captureId = argm.count("c") == 0 ? 0 : std::stoi(argm.at("c"));
    const int nChannels = argm.count("C") == 0 ? 0 : std::stoi(argm.at("C"));

    const int nKeysToCapture = atoi(argv[2]);

    int nKeysHave = 0;
    bool doRecord = true;
    bool doneRecording = false;
    size_t totalSize_bytes = 0;

    TWaveformF waveformF;
    TWaveformI16 waveformI16;
    TKeyPressCollectionT<TSampleI16> keyPresses;

    AudioLogger audioLogger;

    AudioLogger::Callback cbAudio = [&](const auto & frames) {
        for (auto & frame : frames) {
            waveformF.insert(waveformF.end(), frame.begin(), frame.end());
        }

        {
            // apply default filtering, because keypress detection without it is impossible
            auto waveformFFiltered = waveformF;
            ::filter(waveformFFiltered, EAudioFilter::FirstOrderHighPass, kFreqCutoff_Hz, kSampleRate);

            if (convert(waveformFFiltered, waveformI16) == false) {
                printf("Conversion failed\n");
            }
        }

        TWaveformI16 waveformMax;
        TWaveformI16 waveformThreshold;
        if (findKeyPresses(getView(waveformI16, 0), keyPresses, waveformThreshold, waveformMax, 8.0, 512, 2*1024, true) == false) {
            printf("Failed to detect keypresses\n");
        }

        if (nKeysHave < (int) keyPresses.size()) {
            nKeysHave = keyPresses.size();
            printf("Detected %d keys\n", nKeysHave);
        }

        if (nKeysToCapture <= (int) keyPresses.size()) {
            doneRecording = true;
        } else {
            doRecord = true;
        }
    };

    AudioLogger::Parameters parameters;
    parameters.callback = std::move(cbAudio);
    parameters.captureId = captureId;
    parameters.nChannels = nChannels;
    parameters.sampleRate = kSampleRate;
    parameters.filter = EAudioFilter::None;
    //parameters.freqCutoff_Hz = kFreqCutoff_Hz;

    if (audioLogger.install(std::move(parameters)) == false) {
        fprintf(stderr, "Failed to install audio logger\n");
        return -1;
    }

    while (true) {
        if (doRecord) {
            doRecord = false;
            audioLogger.record(0.5f, 0);
        } else {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }

        if (doneRecording) break;
    }

    printf("\n\n Done recording\n");

    std::ofstream fout(argv[1], std::ios::binary);
    if (fout.good() == false) {
        fprintf(stderr, "Failed to open file '%s'\n", argv[1]);
        return -1;
    }

    totalSize_bytes = sizeof(waveformF[0])*waveformF.size();
    fout.write((char *)(waveformF.data()), sizeof(waveformF[0])*waveformF.size());
    fout.close();

    printf("Total data saved: %g MB\n", ((float)(totalSize_bytes)/1024.0f/1024.0f));

    return 0;
}
