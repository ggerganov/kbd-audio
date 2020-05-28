/*! \file guess_qp.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#ifdef __EMSCRIPTEN__
#include "build_timestamp.h"
#include "emscripten/emscripten.h"
#endif

#include "constants.h"
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
#ifdef __EMSCRIPTEN__
    constexpr float kBufferSize_s = 1.0f;
    constexpr uint64_t kSampleRate = 12000;
#else
    constexpr float kBufferSize_s = 0.1f;
    constexpr uint64_t kSampleRate = 48000;
#endif

    constexpr uint64_t kBufferSize_frames = getBufferSize_frames(kSampleRate, kBufferSize_s);

    using ValueCC = float;
    using Offset = int;

    using TKey = int;
    using TKeyWaveform = std::array<AudioLogger::Frame, kBufferSize_frames>;
    using TKeyHistory = std::vector<TKeyWaveform>;

    TKey keyPressed = -1;
    std::map<TKey, TKeyHistory> keySoundHistoryAmpl;
    std::map<TKey, TKeyWaveform> keySoundAverageAmpl;

    int timesToPressQ = 5;
    int timesToPressP = 5;

    bool printStatus = true;
    bool isReadyToPredict = false;

    AudioLogger audioLogger;

    auto calcCC = [](const TKeyWaveform & waveform0, const TKeyWaveform & waveform1, int nSamplesPerFrame, int scmp0, int scmp1, int alignWindow) {
        Offset besto = -1;
        ValueCC bestcc = 0.0f;

        for (int o = -alignWindow; o < alignWindow; ++o) {
            float cc = 0.0f;

            float sum0 = 0.0f, sum02 = 0.0f, sum1 = 0.0f, sum12 = 0.0f, sum01 = 0.0f;
            for (int is = scmp0; is < scmp1; ++is) {
                int is0 = is;
                int f0 = is0/nSamplesPerFrame;
                int s0 = is0 - f0*nSamplesPerFrame;

                int is1 = is + o;
                int f1 = is1/nSamplesPerFrame;
                int s1 = is1 - f1*nSamplesPerFrame;

                auto a0 = waveform0[f0][s0];
                sum0 += a0;
                sum02 += a0*a0;

                auto a1 = waveform1[f1][s1];
                sum1 += a1;
                sum12 += a1*a1;

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

        return std::tuple<ValueCC, Offset>{ bestcc, besto };
    };

    AudioLogger::Callback cbAudio = [&](const AudioLogger::Record & frames) {
        if (frames.size() != keySoundAverageAmpl[keyPressed].size()) {
            printf("Unexpected number of frames - %d, expected - %d. Should never happen\n",
                   (int) frames.size(), (int) keySoundAverageAmpl[keyPressed].size());
            return;
        }

        const int nFrames = frames.size();

        if (isReadyToPredict) {
            TKeyWaveform ampl;
            for (int k = 0; k < nFrames; ++k) {
                ampl[k] = frames[k];
            }

            int nFramesPerWaveform = ampl.size();
            int nSamplesPerFrame = kSamplesPerFrame;
            int centerSample = nFramesPerWaveform*nSamplesPerFrame/2;
            int alignWindow = centerSample/2;

            int scmp0 = centerSample - nSamplesPerFrame/4;
            int scmp1 = centerSample + nSamplesPerFrame/4;

            auto resq = calcCC(keySoundAverageAmpl['q'], ampl, nSamplesPerFrame, scmp0, scmp1, alignWindow);
            auto resp = calcCC(keySoundAverageAmpl['p'], ampl, nSamplesPerFrame, scmp0, scmp1, alignWindow);

            char res = resq > resp ? 'q' : 'p';

            printf("    Prediction: '%c'        ('q' %8.5g vs %8.5g 'p')    (offsets    %5d %5d)\n",
                   res, std::get<0>(resq), std::get<0>(resp), std::get<1>(resq), std::get<1>(resp));

            g_predictedKey = res;
        } else {
            auto & history = keySoundHistoryAmpl[keyPressed];
            history.push_back(TKeyWaveform());
            auto & ampl = history.back();
            for (int k = 0; k < nFrames; ++k) {
                ampl[k] = frames[k];
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
        if (keyPressed == -1) {
            g_predictedKey = -1;
            keyPressed = key;
            audioLogger.record(kBufferSize_s, 3);
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
                int nFramesPerWaveform = history[0].size();
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
                            float acur = std::abs(waveform[iframe][isample]);
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
                    float sum = 0.0f;
                    float sum2 = 0.0f;
                    for (const auto & p : data) {
                        sum += p;
                        sum2 += p*p;
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

                    auto newWaveform = TKeyWaveform();
                    auto & waveform = history[iwaveform];
                    for (int iframe = 0; iframe < nFramesPerWaveform; ++iframe) {
                        for (int isample = 0; isample < nSamplesPerFrame; ++isample) {
                            int icur = iframe*nSamplesPerFrame + isample;
                            int iorg = icur + offset;

                            if (iorg >= 0 && iorg < nFramesPerWaveform*nSamplesPerFrame) {
                                int f = iorg/nSamplesPerFrame;
                                int s = iorg - f*nSamplesPerFrame;
                                newWaveform[iframe][isample] = waveform[f][s];
                            } else {
                                newWaveform[iframe][isample] = 0.0f;
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

                    auto a = waveform0[f][s];
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

                            auto a1 = waveform1[f1][s1];
                            sum1 += a1;
                            sum12 += a1*a1;

                            int is0 = is;
                            int f0 = is0/nSamplesPerFrame;
                            int s0 = is0 - f0*nSamplesPerFrame;

                            auto a0 = waveform0[f0][s0];
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

                    auto newWaveform = TKeyWaveform();
                    for (int iframe = 0; iframe < nFramesPerWaveform; ++iframe) {
                        for (int isample = 0; isample < nSamplesPerFrame; ++isample) {
                            int icur = iframe*nSamplesPerFrame + isample;
                            int iorg = icur + besto;

                            if (iorg >= 0 && iorg < nFramesPerWaveform*nSamplesPerFrame) {
                                int f = iorg/nSamplesPerFrame;
                                int s = iorg - f*nSamplesPerFrame;
                                newWaveform[iframe][isample] = waveform1[f][s];
                            } else {
                                newWaveform[iframe][isample] = 0.0f;
                            }
                        }
                    }

                    waveform1 = std::move(newWaveform);
                }

                printf("    - Calculating average waveform\n");
                auto & avgWaveform = keySoundAverageAmpl[key];
                for (auto & f : avgWaveform) f.fill(0.0f);
                for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                    auto & waveform = history[iwaveform];
                    for (int iframe = 0; iframe < nFramesPerWaveform; ++iframe) {
                        for (int isample = 0; isample < nSamplesPerFrame; ++isample) {
                            avgWaveform[iframe][isample] += waveform[iframe][isample];
                        }
                    }
                }

                {
                    float norm = 1.0f/(nFramesPerWaveform*nSamplesPerFrame);
                    for (int iframe = 0; iframe < nFramesPerWaveform; ++iframe) {
                        for (int isample = 0; isample < nSamplesPerFrame; ++isample) {
                            avgWaveform[iframe][isample] *= norm;
                        }
                    }
                }

                printf("\n");
            };

            trainKey('q');
            trainKey('p');
            isReadyToPredict = true;

            printf("[+] Ready to predict. Keep pressing 'q' or 'p' and the program will guess which key was pressed\n");
            printf("    based on the captured audio from the microphone.\n");
            printf("[+] Predicting\n");
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
