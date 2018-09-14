/*! \file guess_qp.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#ifdef __EMSCRIPTEN__
#include "emscripten/emscripten.h"
#endif

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

std::function<void()> g_update;
std::function<void(int)> g_handleKey;

extern "C" {
    void keyPressedCallback(int key) {
        g_handleKey(key);
    }
}

void update() {
    g_update();
}

int main(int, char**) {
    using ValueCC = float;
    using Offset = int;

    using TKey = int;
    using TKeyWaveform = std::array<AudioLogger::Frame, 2*AudioLogger::kBufferSize_frames - 1>;
    using TKeyHistory = std::vector<TKeyWaveform>;

    TKey keyPressed = -1;
    std::map<TKey, TKeyHistory> keySoundHistoryAmpl;
    std::map<TKey, TKeyWaveform> keySoundAverageAmpl;

    int timesToPressQ = 5;
    int timesToPressP = 5;

    bool printStatus = true;
    bool isReadyToPredict = false;

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

    AudioLogger audioLogger;
    AudioLogger::Callback cbAudio = [&](const AudioLogger::Record & frames) {
        if (frames.size() != keySoundAverageAmpl[keyPressed].size()) {
            printf("Unexpected number of frames - %d. Should never happen\n", (int) frames.size());
            return;
        }

        const int nFrames = frames.size();

        if (isReadyToPredict) {
            TKeyWaveform ampl;
            for (int k = 0; k < nFrames; ++k) {
                ampl[k] = frames[k];
            }

            int nFramesPerStroke = ampl.size();
            int nSamplesPerFrame = AudioLogger::kSamplesPerFrame;
            int centerSample = nFramesPerStroke*nSamplesPerFrame/2;
            int alignWindow = nSamplesPerFrame;

            int scmp0 = centerSample - alignWindow;
            int scmp1 = centerSample + alignWindow;

            auto resq = calcCC(keySoundAverageAmpl['q'], ampl, nSamplesPerFrame, scmp0, scmp1, 2*alignWindow);
            auto resp = calcCC(keySoundAverageAmpl['p'], ampl, nSamplesPerFrame, scmp0, scmp1, 2*alignWindow);

            char res = resq > resp ? 'q' : 'p';

            printf("    Prediction: '%c'        ('q' %8.5g vs %8.5g 'p')    (offsets    %5d %5d)\n",
                   res, std::get<0>(resq), std::get<0>(resp), std::get<1>(resq), std::get<1>(resp));
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
            }
            if (keyPressed == 'p' && timesToPressP > 0) {
                --timesToPressP;
                printStatus = true;
            }
        }

        keyPressed = -1;
    };

    if (audioLogger.install(cbAudio) == false) {
        fprintf(stderr, "Failed to install audio logger\n");
        return -1;
    }

    g_handleKey = [&](int key) {
        if (keyPressed == -1) {
            keyPressed = key;
            audioLogger.record();
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
            std::this_thread::sleep_for(std::chrono::milliseconds(1));

            return;
        }

        if (isReadyToPredict == false) {
            printf("[+] Training\n");

            auto trainKey = [&](TKey key) {
                auto & history = keySoundHistoryAmpl[key];

                int nStrokes = history.size();
                int nFramesPerStroke = history[0].size();
                int nSamplesPerFrame = AudioLogger::kSamplesPerFrame;

                printf("    - Training key '%c'\n", key);
                printf("    - History size = %d key strokes\n", nStrokes);
                printf("    - Frames per key stroke     = %d\n", nFramesPerStroke);
                printf("    - Total frames available    = %d\n", nStrokes*nFramesPerStroke);
                printf("    - Samples per frame         = %d\n", nSamplesPerFrame);
                printf("    - Total samples available   = %d\n", nStrokes*nFramesPerStroke*nSamplesPerFrame);

                printf("    - Estimating stroke peaks ...\n");
                std::vector<int> peakSum;
                std::vector<int> peakMax;

                peakSum.clear();
                peakMax.clear();
                for (int istroke = 0; istroke < nStrokes; ++istroke) {
                    int isum = -1;
                    float asum = 0.0f;
                    float aisum = 0.0f;

                    int imax = -1;
                    float amax = 0.0f;

                    const auto & stroke = history[istroke];

                    for (int iframe = 0; iframe < nFramesPerStroke; ++iframe) {
                        for (int isample = 0; isample < nSamplesPerFrame; ++isample) {
                            int icur = iframe*nSamplesPerFrame + isample;
                            float acur = std::abs(stroke[iframe][isample]);
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

                int centerSample = nFramesPerStroke*nSamplesPerFrame/2;

                printf("    - Centering strokes at sample %d\n", centerSample);
                for (int istroke = 0; istroke < nStrokes; ++istroke) {
                    int offset = peakUsed[istroke] - centerSample;
                    printf("        Offset for stroke %d = %d\n", istroke, offset);

                    auto newStroke = TKeyWaveform();
                    auto & stroke = history[istroke];
                    for (int iframe = 0; iframe < nFramesPerStroke; ++iframe) {
                        for (int isample = 0; isample < nSamplesPerFrame; ++isample) {
                            int icur = iframe*nSamplesPerFrame + isample;
                            int iorg = icur + offset;

                            if (iorg >= 0 && iorg < nFramesPerStroke*nSamplesPerFrame) {
                                int f = iorg/nSamplesPerFrame;
                                int s = iorg - f*nSamplesPerFrame;
                                newStroke[iframe][isample] = stroke[f][s];
                            } else {
                                newStroke[iframe][isample] = 0.0f;
                            }
                        }
                    }

                    stroke = std::move(newStroke);
                }

                int alignToStroke = 0;
                int alignWindow = nSamplesPerFrame;
                printf("    - Aligning all strokes to stroke %d using cross correlation\n", alignToStroke);
                printf("      Align window = %d\n", alignWindow);

                int scmp0 = centerSample - alignWindow;
                int scmp1 = centerSample + alignWindow;

                float sum0 = 0.0f;
                float sum02 = 0.0f;

                const auto & stroke0 = history[alignToStroke];

                for (int is = scmp0; is < scmp1; ++is) {
                    int f = is/nSamplesPerFrame;
                    int s = is - f*nSamplesPerFrame;

                    auto a = stroke0[f][s];
                    sum0 += a;
                    sum02 += a*a;
                }

                for (int istroke = 0; istroke < nStrokes; ++istroke) {
                    if (istroke == alignToStroke) continue;

                    auto & stroke1 = history[istroke];

                    int besto = 0;
                    float bestcc = 0.0f;

                    for (int o = -alignWindow; o < alignWindow; ++o) {
                        float cc = 0.0f;

                        float sum1 = 0.0f, sum12 = 0.0f, sum01 = 0.0f;
                        for (int is = scmp0; is < scmp1; ++is) {
                            int is0 = is;
                            int f0 = is0/nSamplesPerFrame;
                            int s0 = is0 - f0*nSamplesPerFrame;

                            int is1 = is + o;
                            int f1 = is1/nSamplesPerFrame;
                            int s1 = is1 - f1*nSamplesPerFrame;

                            auto a1 = stroke1[f1][s1];
                            sum1 += a1;
                            sum12 += a1*a1;

                            auto a0 = stroke0[f0][s0];
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

                    printf("        Best offset for stroke %d = %d (cc = %g)\n", istroke, besto, bestcc);

                    auto newStroke = TKeyWaveform();
                    auto & stroke = history[istroke];
                    for (int iframe = 0; iframe < nFramesPerStroke; ++iframe) {
                        for (int isample = 0; isample < nSamplesPerFrame; ++isample) {
                            int icur = iframe*nSamplesPerFrame + isample;
                            int iorg = icur + besto;

                            if (iorg >= 0 && iorg < nFramesPerStroke*nSamplesPerFrame) {
                                int f = iorg/nSamplesPerFrame;
                                int s = iorg - f*nSamplesPerFrame;
                                newStroke[iframe][isample] = stroke[f][s];
                            } else {
                                newStroke[iframe][isample] = 0.0f;
                            }
                        }
                    }

                    stroke = std::move(newStroke);
                }

                printf("    - Calculating average waveform\n");
                auto & avgWaveform = keySoundAverageAmpl[key];
                for (auto & f : avgWaveform) f.fill(0.0f);
                for (int istroke = 0; istroke < nStrokes; ++istroke) {
                    auto & stroke = history[istroke];
                    for (int iframe = 0; iframe < nFramesPerStroke; ++iframe) {
                        for (int isample = 0; isample < nSamplesPerFrame; ++isample) {
                            avgWaveform[iframe][isample] += stroke[iframe][isample];
                        }
                    }
                }

                {
                    float norm = 1.0f/(nFramesPerStroke*nSamplesPerFrame);
                    for (int iframe = 0; iframe < nFramesPerStroke; ++iframe) {
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

            printf("[+] Predicting\n");
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    };

    printf("[+] Collecting training data\n");
#ifdef __EMSCRIPTEN__
    emscripten_set_main_loop(update, 0, 1);
#else
    while (true) {
        update();
    }
#endif

    return 0;
}
