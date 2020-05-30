/*! \file guess_qp2.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#ifdef __EMSCRIPTEN__
#include "build_timestamp.h"
#include "emscripten/emscripten.h"
#endif

#include "constants.h"
#include "common.h"
#include "audio_logger.h"

#include <stdio.h>
#include <termios.h>
#include <unistd.h>

#include <map>
#include <cmath>
#include <string>
#include <chrono>
#include <thread>
#include <vector>

static int g_predictedKey = -1;
static bool g_isInitialized = false;

static std::function<int()> g_init;
static std::function<void()> g_update;
static std::function<void(int)> g_handleKey;

int init() {
    if (g_isInitialized) return 1;

    return g_init();
}

void update() {
    if (g_isInitialized == false) return;

    g_update();
}

using TKeyWaveform = TKeyWaveformF;
using TKeyHistory = TKeyHistoryF;

// JS interface
extern "C" {
    int doInit() {
        return init();
    }

    void keyPressedCallback(int key) {
        g_handleKey(key);
    }

    int getPredictedKey() {
        return g_predictedKey;
    }
}

int main(int, char**) {
    TKey keyPressed = -1;
    std::map<TKey, TKeyHistory> keySoundHistoryAmpl;
    std::map<TKey, TKeyWaveform> keySoundAverageAmpl;

    int timesToPressQ = 5;
    int timesToPressP = 5;

    bool doRecord = false;
    bool printStatus = true;
    bool isReadyToPredict = false;

    // rig buffer
    int rbBegin = 0;
    float rbAverage = 0.0f;
    std::array<float, kBkgrRingBufferSize> rbSamples;
    rbSamples.fill(0.0f);

    int lastkey = -1;
    double lastcc = -1.0f;
    float thresholdCC = 0.35f;

    AudioLogger audioLogger;

    AudioLogger::Callback cbAudio = [&](const AudioLogger::Record & frames) {
        const int nFrames = frames.size();

        if (isReadyToPredict) {

            std::vector<int> positionsToPredict;

            {
                float amax = 0.0f;
                for (int f = 0; f < (int) frames.size(); ++f) {
                    for (int s = 0; s < (int) frames[f].size(); s += kBkgrStep_samples) {
                        rbAverage *= rbSamples.size();
                        rbAverage -= rbSamples[rbBegin];
                        auto acur = std::abs(frames[f][s]);
                        rbSamples[rbBegin] = acur;
                        if (acur > amax) amax = acur;
                        rbAverage += acur;
                        rbAverage /= rbSamples.size();
                        if (++rbBegin >= (int) rbSamples.size()) rbBegin = 0;
                    }
                }

                int skip_samples = 0;
                for (int f = nFrames/2 - 1; f <= nFrames/2 + 1; ++f) {
                    for (int s = 0; s < (int) frames[f].size(); ++s) {
                        if (s + skip_samples >= (int) frames[f].size()) {
                            skip_samples -= frames[f].size() - s;
                            s += skip_samples;
                            continue;
                        } else {
                            s += skip_samples;
                            skip_samples = 0;
                        }
                        auto acur = frames[f][s];
                        if (acur > 10.0f*rbAverage) {
                            skip_samples = kKeyDuration_samples;
                            positionsToPredict.push_back(f*kSamplesPerFrame + s);
                            //printf("Key press detected\n");
                        }
                    }
                }

                //printf("Average = %10.8f, max = %10.8f\n", rbAverage, amax);
            }

            if (positionsToPredict.size() > 0) {
                TKeyWaveform ampl(nFrames*kSamplesPerFrame);
                for (int k = 0; k < nFrames; ++k) {
                    std::copy(frames[k].begin(), frames[k].end(), ampl.begin() + k*kSamplesPerFrame);
                }

                int alignWindow = 0.10*kSamplesPerWaveformTrain;

                for (int ipos = 0; ipos < (int) positionsToPredict.size() ; ++ipos) {
                    int scmp0 = positionsToPredict[ipos] - 0.25*kSamplesPerWaveformTrain;
                    int scmp1 = positionsToPredict[ipos] + 0.25*kSamplesPerWaveformTrain;

                    auto resq = findBestCC(keySoundAverageAmpl['q'], ampl, scmp0, scmp1, alignWindow);
                    auto resp = findBestCC(keySoundAverageAmpl['p'], ampl, scmp0, scmp1, alignWindow);

                    char res = std::get<0>(resq) > std::get<0>(resp) ? 'q' : 'p';

                    auto curcc = std::max(std::get<0>(resp), std::get<0>(resq));
                    if (curcc > thresholdCC) {
                        if (lastkey != res || lastcc != curcc) {
                            printf("    Prediction: '%c'        ('q' %8.5g vs %8.5g 'p')    (offsets    %5ld %5ld)\n",
                                   res, std::get<0>(resq), std::get<0>(resp), std::get<1>(resq), std::get<1>(resp));
                            lastkey = res;
                            lastcc = curcc;
                        }
                    }

                    g_predictedKey = res;
                }
            }

            doRecord = true;
        } else {
            auto & history = keySoundHistoryAmpl[keyPressed];
            history.push_back(TKeyWaveform(nFrames*kSamplesPerFrame));
            auto & ampl = history.back();
            for (int k = 0; k < nFrames; ++k) {
                std::copy(frames[k].begin(), frames[k].end(), ampl.begin() + k*kSamplesPerFrame);
            }

            if (keyPressed == 'q' && timesToPressQ > 0) {
                --timesToPressQ;
                printStatus = true;
            } else if (keyPressed == 'p' && timesToPressP > 0) {
                --timesToPressP;
                printStatus = true;
            }
        }

        keyPressed = -1;
    };

    g_init = [&]() {
        AudioLogger::Parameters parameters;
        parameters.callback = std::move(cbAudio);
        parameters.captureId = 0;
        parameters.nChannels = 1;
        parameters.sampleRate = kSampleRate;
        parameters.freqCutoff_Hz = kFreqCutoff_Hz;

        if (audioLogger.install(std::move(parameters)) == false) {
            fprintf(stderr, "Failed to install audio logger\n");
            return -1;
        }

        printf("[+] Collecting training data\n");
        g_isInitialized = true;
        return 0;
    };

    g_handleKey = [&](int key) {
        if (keyPressed == -1 && isReadyToPredict == false) {
            g_predictedKey = -1;
            keyPressed = key;
            audioLogger.record(kBufferSizeTrain_s, 3);
        }
    };

#ifdef __EMSCRIPTEN__
#else
    std::thread keyReader = std::thread([&]() {
        struct termios oldt, newt;
        tcgetattr ( STDIN_FILENO, &oldt );
        newt = oldt;
        newt.c_lflag &= ~( ICANON | ECHO );
        tcsetattr ( STDIN_FILENO, TCSANOW, &newt );
        while (true) {
            int key = getchar();
            g_handleKey(key);
        }
        tcsetattr ( STDIN_FILENO, TCSANOW, &oldt );
    });
#endif

    g_update = [&]() {
        if (timesToPressQ > 0 || timesToPressP > 0) {
            if (printStatus) {
                if (timesToPressQ > 0) {
                    printf("    - press the letter 'q' %d more times\n", timesToPressQ);
                } else if (timesToPressP > 0) {
                    printf("    - press the letter 'p' %d more times\n", timesToPressP);
                }
                printStatus = false;
            }

            return;
        }

        if (isReadyToPredict == false) {
            printf("[+] Training\n");

            auto trainKey = [&](TKey key) {
                auto & history = keySoundHistoryAmpl[key];

                int nWaveforms = history.size();
                int nFramesPerWaveform = history[0].size()/kSamplesPerFrame;
                int nSamplesPerFrame = kSamplesPerFrame;

                printf("    - Training key '%c'\n", key);
                printf("    - History size = %d key waveforms\n", nWaveforms);
                printf("    - Frames per key waveform   = %d\n", nFramesPerWaveform);
                printf("    - Total frames available    = %d\n", nWaveforms*nFramesPerWaveform);
                printf("    - Samples per frame         = %d\n", nSamplesPerFrame);
                printf("    - Total samples available   = %d\n", nWaveforms*nFramesPerWaveform*nSamplesPerFrame);

                printf("    - Estimating waveform peaks ...\n");
                std::vector<int> peakSum;
                std::vector<int> peakMax;

                peakSum.clear();
                peakMax.clear();
                for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                    int isum = -1;
                    float asum = 0.0f;
                    float aisum = 0.0f;

                    int imax = -1;
                    float amax = 0.0f;

                    const auto & waveform = history[iwaveform];

                    for (int iframe = 0; iframe < nFramesPerWaveform; ++iframe) {
                        for (int isample = 0; isample < nSamplesPerFrame; ++isample) {
                            int icur = iframe*nSamplesPerFrame + isample;
                            float acur = std::abs(waveform[iframe*kSamplesPerFrame + isample]);
                            float acur2 = acur*acur;

                            asum += acur2;
                            aisum += acur2*icur;

                            if (acur > amax) {
                                amax = acur;
                                imax = icur;
                            }
                        }
                    }

                    isum = aisum/asum;

                    peakSum.push_back(isum);
                    peakMax.push_back(imax);
                    printf("        Estimated peak: %d (method - sum), %d (method - max)\n", isum, imax);
                }

                auto calcStdev = [](const std::vector<int> & data) {
                    double sum = 0.0f;
                    double sum2 = 0.0f;
                    for (const auto & p : data) {
                        int64_t v = p;
                        sum += v;
                        sum2 += v*v;
                    }
                    sum /= data.size();
                    sum2 /= data.size();
                    return sqrt(sum2 - sum*sum);
                };

                float stdevSum = calcStdev(peakSum);
                float stdevMax = calcStdev(peakMax);

                printf("    - Stdev of estimated peaks: %g (sum) vs %g (max)\n", stdevSum, stdevMax);

                const auto & peakUsed = peakMax;
                printf("    - Using 'max' estimation\n");

                int centerSample = nFramesPerWaveform*nSamplesPerFrame/2;

                printf("    - Centering waveforms at sample %d\n", centerSample);
                for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                    int offset = peakUsed[iwaveform] - centerSample;
                    printf("        Offset for waveform %d = %d\n", iwaveform, offset);

                    auto & waveform = history[iwaveform];
                    auto newWaveform = TKeyWaveform(waveform.size());
                    for (int iframe = 0; iframe < nFramesPerWaveform; ++iframe) {
                        for (int isample = 0; isample < nSamplesPerFrame; ++isample) {
                            int icur = iframe*nSamplesPerFrame + isample;
                            int iorg = icur + offset;

                            if (iorg >= 0 && iorg < nFramesPerWaveform*nSamplesPerFrame) {
                                int f = iorg/nSamplesPerFrame;
                                int s = iorg - f*nSamplesPerFrame;
                                newWaveform[iframe*kSamplesPerFrame + isample] = waveform[f*kSamplesPerFrame + s];
                            } else {
                                newWaveform[iframe*kSamplesPerFrame + isample] = 0.0f;
                            }
                        }
                    }

                    waveform = std::move(newWaveform);
                }

                int alignToWaveform = nWaveforms/2;
                int alignWindow = centerSample/2;
                printf("    - Aligning all waveforms to waveform %d using cross correlation\n", alignToWaveform);
                printf("      Align window = %d\n", alignWindow);

                int scmp0 = centerSample - nSamplesPerFrame/4;
                int scmp1 = centerSample + nSamplesPerFrame/4;

                float sum0 = 0.0f;
                float sum02 = 0.0f;

                const auto & waveform0 = history[alignToWaveform];

                for (int is = scmp0; is < scmp1; ++is) {
                    int f = is/nSamplesPerFrame;
                    int s = is - f*nSamplesPerFrame;

                    auto a = waveform0[f*kSamplesPerFrame + s];
                    sum0 += a;
                    sum02 += a*a;
                }

                for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                    if (iwaveform == alignToWaveform) continue;

                    auto & waveform1 = history[iwaveform];

                    int besto = 0;
                    float bestcc = 0.0f;

                    for (int o = -alignWindow; o < alignWindow; ++o) {
                        float cc = 0.0f;

                        float sum1 = 0.0f, sum12 = 0.0f, sum01 = 0.0f;
                        for (int is = scmp0; is < scmp1; ++is) {
                            int is1 = is + o;
                            int f1 = is1/nSamplesPerFrame;
                            int s1 = is1 - f1*nSamplesPerFrame;

                            auto a1 = waveform1[f1*kSamplesPerFrame + s1];
                            sum1 += a1;
                            sum12 += a1*a1;

                            int is0 = is;
                            int f0 = is0/nSamplesPerFrame;
                            int s0 = is0 - f0*nSamplesPerFrame;

                            auto a0 = waveform0[f0*kSamplesPerFrame + s0];
                            sum01 += a0*a1;
                        }

                        int ncc = scmp1 - scmp0;
                        {
                            float nom = sum01*ncc - sum0*sum1;
                            float den2a = sum02*ncc - sum0*sum0;
                            float den2b = sum12*ncc - sum1*sum1;
                            cc = (nom)/(sqrt(den2a*den2b));
                        }

                        if (cc > bestcc) {
                            besto = o;
                            bestcc = cc;
                        }
                    }

                    printf("        Best offset for waveform %d = %d (cc = %g)\n", iwaveform, besto, bestcc);

                    auto newWaveform = TKeyWaveform(nFramesPerWaveform*nSamplesPerFrame);
                    for (int iframe = 0; iframe < nFramesPerWaveform; ++iframe) {
                        for (int isample = 0; isample < nSamplesPerFrame; ++isample) {
                            int icur = iframe*nSamplesPerFrame + isample;
                            int iorg = icur + besto;

                            if (iorg >= 0 && iorg < nFramesPerWaveform*nSamplesPerFrame) {
                                int f = iorg/nSamplesPerFrame;
                                int s = iorg - f*nSamplesPerFrame;
                                newWaveform[iframe*kSamplesPerFrame + isample] = waveform1[f*kSamplesPerFrame + s];
                            } else {
                                newWaveform[iframe*kSamplesPerFrame + isample] = 0.0f;
                            }
                        }
                    }

                    waveform1 = std::move(newWaveform);
                }

                printf("    - Calculating average waveform\n");
                auto & avgWaveform = keySoundAverageAmpl[key];
                avgWaveform.resize(nFramesPerWaveform*nSamplesPerFrame);
                for (auto & f : avgWaveform) f = 0.0f;
                for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                    auto & waveform = history[iwaveform];
                    for (int iframe = 0; iframe < nFramesPerWaveform; ++iframe) {
                        for (int isample = 0; isample < nSamplesPerFrame; ++isample) {
                            avgWaveform[iframe*kSamplesPerFrame + isample] += waveform[iframe*kSamplesPerFrame + isample];
                        }
                    }
                }

                {
                    float norm = 1.0f/(nFramesPerWaveform*nSamplesPerFrame);
                    for (int iframe = 0; iframe < nFramesPerWaveform; ++iframe) {
                        for (int isample = 0; isample < nSamplesPerFrame; ++isample) {
                            avgWaveform[iframe*kSamplesPerFrame + isample] *= norm;
                        }
                    }
                }

                printf("\n");
            };

            trainKey('q');
            trainKey('p');
            isReadyToPredict = true;
            doRecord = true;

            printf("[+] Ready to predict. Keep pressing 'q' or 'p' and the program will guess which key was pressed\n");
            printf("    based on the captured audio from the microphone.\n");
            printf("[+] Predicting\n");
        }

        if (doRecord) {
            doRecord = false;
            audioLogger.record(kBufferSizeTrain_s, getBufferSize_frames(kSampleRate, kBufferSizeTrain_s) - 1);
        }
    };

#ifdef __EMSCRIPTEN__
    printf("Build time: %s\n", BUILD_TIMESTAMP);
    printf("Press the Init button to start\n");
    emscripten_set_main_loop(update, 60, 1);
#else
    init();
    while (true) {
        update();
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
#endif

    return 0;
}
