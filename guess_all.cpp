/*! \file guess_all.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "audio_logger.h"

#include <map>
#include <cmath>
#include <string>
#include <chrono>
#include <thread>
#include <vector>
#include <deque>
#include <fstream>

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

int main(int argc, char ** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.kbd [input2.kbd ...]\n", argv[0]);
        return -127;
    }

    constexpr float kBufferSize_s = 0.075f;
    constexpr uint64_t kSampleRate = 96000;

    constexpr uint64_t kRingBufferSize = 16*1024;
    constexpr int bkgrStep_samples = 7;
    constexpr int keyDuration_samples = 0.050f*kSampleRate;

    constexpr uint64_t kBufferSize_frames = 2*AudioLogger::getBufferSize_frames(kSampleRate, kBufferSize_s) - 1;

    std::map<int, std::ifstream> fins;
    for (int i = 0; i < argc - 1; ++i) {
        printf("Opening file '%s'\n", argv[i + 1]);
        fins[i] = std::ifstream(argv[i + 1], std::ios::binary);

        {
            int bufferSize_frames = 1;
            fins[i].read((char *)(&bufferSize_frames), sizeof(bufferSize_frames));
            if (bufferSize_frames != kBufferSize_frames) {
                printf("Buffer size in file (%d) does not match the expected one (%d)\n", bufferSize_frames, (int) kBufferSize_frames);
                return -1;
            }
        }
    }

    using ValueCC = float;
    using Offset = int;

    using TKey = int;
    using TKeyWaveform = std::array<AudioLogger::Frame, kBufferSize_frames>;
    using TKeyHistory = std::vector<TKeyWaveform>;

    TKey keyPressed = -1;
    std::map<TKey, TKeyHistory> keySoundHistoryAmpl;
    std::map<TKey, TKeyWaveform> keySoundAverageAmpl;

    bool doRecord = false;
    bool printStatus = true;
    bool isReadyToPredict = false;
    bool processingInput = true;

    int curFile = 0;

    // rig buffer
    int rbBegin = 0;
    float rbAverage = 0.0f;
    std::array<float, kRingBufferSize> rbSamples;
    rbSamples.fill(0.0f);

    AudioLogger audioLogger;

    auto calcCC = [](const TKeyWaveform & waveform0, const TKeyWaveform & waveform1, int nSamplesPerFrame, int scmp0, int scmp1, int alignWindow) {
        int is00 = waveform0.size()*nSamplesPerFrame/2 - (scmp1 - scmp0)/2;

        Offset besto = -1;
        ValueCC bestcc = -1.0f;

        float sum0 = 0.0f, sum02 = 0.0f;
        for (int is = 0; is < scmp1 - scmp0; ++is) {
            int is0 = is00 + is;
            int f0 = is0/nSamplesPerFrame;
            int s0 = is0 - f0*nSamplesPerFrame;

            auto a0 = waveform0[f0][s0];
            sum0 += a0;
            sum02 += a0*a0;
        }

        int nWorkers = 4;
        std::mutex mutex;
        std::vector<std::thread> workers(nWorkers);
        for (int i = 0; i < workers.size(); ++i) {
            auto & worker = workers[i];
            worker = std::thread([&, i]() {
                Offset cbesto = -1;
                ValueCC cbestcc = -1.0f;

                for (int o = -alignWindow + i; o < alignWindow; o += nWorkers) {
                    float cc = -1.0f;

                    float sum1 = 0.0f, sum12 = 0.0f, sum01 = 0.0f;
                    for (int is = 0; is < scmp1 - scmp0; ++is) {
                        int is0 = is00 + is;
                        int f0 = is0/nSamplesPerFrame;
                        int s0 = is0 - f0*nSamplesPerFrame;

                        int is1 = is + scmp0 + o;
                        int f1 = is1/nSamplesPerFrame;
                        int s1 = is1 - f1*nSamplesPerFrame;

                        auto a0 = waveform0[f0][s0];
                        auto a1 = waveform1[f1][s1];
                        sum1 += a1;
                        sum12 += a1*a1;

                        sum01 += a0*a1;
                    }

                    int ncc = (scmp1 - scmp0);
                    {
                        float nom = sum01*ncc - sum0*sum1;
                        float den2a = sum02*ncc - sum0*sum0;
                        float den2b = sum12*ncc - sum1*sum1;
                        cc = (nom)/(sqrt(den2a*den2b));
                    }

                    if (cc > cbestcc) {
                        cbesto = o;
                        cbestcc = cc;
                    }
                }
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    if (cbestcc > bestcc) {
                        bestcc = cbestcc;
                        besto = cbesto;
                    }
                }
            });
        }
        for (auto & worker : workers) worker.join();
        return std::tuple<ValueCC, Offset>{ bestcc, besto };
    };

    struct WorkData {
        TKeyWaveform ampl;
        std::vector<int> positionsToPredict;
    };

    std::mutex mutex;
    std::deque<WorkData> workQueue;
    std::thread worker([&]() {
        int lastkey = -1;
        float lastcc = -1.0f;

        while(true) {
            bool process = false;
            WorkData workData;
            {
                std::lock_guard<std::mutex> lock(mutex);
                while (workQueue.size() > 30) {
                    workQueue.pop_front();
                    printf("pop\n");
                }

                if (workQueue.size() > 0) {
                    workData = std::move(workQueue.front());
                    workQueue.pop_front();
                    process = true;
                }
            }
            if (process) {
                int nSamplesPerFrame = AudioLogger::kSamplesPerFrame;
                const auto & ampl = workData.ampl;
                const auto & positionsToPredict = workData.positionsToPredict;

                int nFramesPerWaveform = ampl.size();
                int alignWindow = 2*nSamplesPerFrame;

                for (int ipos = 0; ipos < positionsToPredict.size() ; ++ipos) {
                    int scmp0 = positionsToPredict[ipos] - nSamplesPerFrame/2;
                    int scmp1 = positionsToPredict[ipos] + nSamplesPerFrame/2;

                    char res = -1;
                    float maxcc = -1.0f;
                    for (const auto & ka : keySoundAverageAmpl) {
                        auto curcc = calcCC(keySoundAverageAmpl[ka.first], ampl, nSamplesPerFrame, scmp0, scmp1, alignWindow);
                        if (std::get<0>(curcc) > maxcc) {
                            res = ka.first;
                            maxcc = std::get<0>(curcc);
                        }
                    }

                    if (maxcc > 0.50f) {
                        if (lastkey != res || lastcc != maxcc) {
                            printf("    Prediction: '%c'        (%8.5g)\n", res, maxcc);
                        }
                        lastkey = res;
                        lastcc = maxcc;
                    }

                    g_predictedKey = res;
                }
            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
    });

    AudioLogger::Callback cbAudio = [&](const AudioLogger::Record & frames) {
        if (frames.size() != keySoundAverageAmpl[keyPressed].size()) {
            printf("Unexpected number of frames - %d, expected - %d. Should never happen\n",
                   (int) frames.size(), (int) keySoundAverageAmpl[keyPressed].size());
            return;
        }

        const int nFrames = frames.size();

        if (isReadyToPredict) {

            int nSamplesPerFrame = AudioLogger::kSamplesPerFrame;

            std::vector<int> positionsToPredict;

            {
                float amax = 0.0f;
                for (int f = 0; f < frames.size(); ++f) {
                    for (int s = 0; s < frames[f].size(); s += bkgrStep_samples) {
                        rbAverage *= rbSamples.size();
                        rbAverage -= rbSamples[rbBegin];
                        auto acur = std::abs(frames[f][s]);
                        rbSamples[rbBegin] = acur;
                        if (acur > amax) amax = acur;
                        rbAverage += acur;
                        rbAverage /= rbSamples.size();
                        if (++rbBegin >= rbSamples.size()) rbBegin = 0;
                    }
                }

                int skip_samples = 0;
                int nFrames = frames.size();
                int nFrames2 = std::max(1, nFrames/2);
                for (int f = nFrames2 - nFrames2/2; f <= nFrames2 + nFrames2/2; ++f) {
                    for (int s = 0; s < frames[f].size(); ++s) {
                        if (s + skip_samples >= frames[f].size()) {
                            skip_samples -= frames[f].size() - s;
                            s += skip_samples;
                            continue;
                        } else {
                            s += skip_samples;
                            skip_samples = 0;
                        }
                        auto acur = frames[f][s];
                        if (acur > 10.0f*rbAverage) {
                            skip_samples = keyDuration_samples;
                            positionsToPredict.push_back(f*nSamplesPerFrame + s);
                            //printf("Key press detected\n");
                        }
                    }
                }

                //printf("Average = %10.8f, max = %10.8f\n", rbAverage, amax);
            }

            if (positionsToPredict.size() > 0) {
                WorkData workData;
                for (int k = 0; k < nFrames; ++k) {
                    workData.ampl[k] = frames[k];
                }
                workData.positionsToPredict = positionsToPredict;

                {
                    std::lock_guard<std::mutex> lock(mutex);
                    workQueue.push_back(std::move(workData));
                }
            }

            doRecord = true;
        } else {
            auto & history = keySoundHistoryAmpl[keyPressed];
            history.push_back(TKeyWaveform());
            auto & ampl = history.back();
            for (int k = 0; k < nFrames; ++k) {
                ampl[k] = frames[k];
            }
        }

        keyPressed = -1;
    };

    g_init = [&]() {
        if (audioLogger.install(kSampleRate, cbAudio) == false) {
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
            audioLogger.record(kBufferSize_s);
        }
    };

    g_update = [&]() {
        if (processingInput) {
            if (keyPressed == -1) {
                AudioLogger::Frame frame;
                AudioLogger::Record record;
                fins[curFile].read((char *)(&keyPressed), sizeof(keyPressed));
                if (fins[curFile].eof()) {
                    ++curFile;
                    if (curFile >= fins.size()) {
                        processingInput = false;
                    }
                } else {
                    printf("%c", keyPressed);
                    fflush(stdout);
                    for (int i = 0; i < kBufferSize_frames; ++i) {
                        fins[curFile].read((char *)(frame.data()), sizeof(AudioLogger::Sample)*frame.size());
                        record.push_back(frame);
                    }
                    cbAudio(record);
                }
            }
            return;
        }

        if (isReadyToPredict == false) {
            printf("[+] Training\n");

            std::vector<TKey> failedToTrain;

            auto trainKey = [&](TKey key) {
                auto & history = keySoundHistoryAmpl[key];

                int nWaveforms = history.size();
                int nFramesPerWaveform = history[0].size();
                int nSamplesPerFrame = AudioLogger::kSamplesPerFrame;

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
                    printf("        Offset for waveform %-4d = %-4d\n", iwaveform, offset);

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

                bool aligned = false;
                std::map<int, float> aligncc;
                for (int alignToWaveform = 0; alignToWaveform < nWaveforms; ++alignToWaveform) {
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

                    int nconf = 0;
                    float sumcc = 0.0f;

                    auto curHistory = history;
                    for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                        aligncc[iwaveform] = 1.0f;
                        if (iwaveform == alignToWaveform) continue;

                        auto & waveform1 = curHistory[iwaveform];

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

                        if (bestcc > 0.65f) ++nconf;
                        sumcc += bestcc;
                        aligncc[iwaveform] = bestcc;
                        printf("        Best offset for waveform %-4d = %-4d (cc = %g)\n", iwaveform, besto, bestcc);

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

                    if (nconf > nWaveforms/2) {
                        history = curHistory;
                        aligned = true;
                        break;
                    }
                }

                if (aligned == false) {
                    failedToTrain.push_back(key);
                }

                printf("    - Calculating average waveform\n");
                auto & avgWaveform = keySoundAverageAmpl[key];
                for (auto & f : avgWaveform) f.fill(0.0f);
                for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                    if (aligncc[iwaveform] < 0.65f) continue;
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

            for (const auto & kh : keySoundHistoryAmpl) {
                if (kh.second.size() > 2) {
                    trainKey(kh.first);
                } else {
                    failedToTrain.push_back(kh.first);
                }
            }
            printf("Failed to train the following keys: ");
            for (auto & k : failedToTrain) printf("'%c' ", k);
            printf("\n");
            isReadyToPredict = true;
            doRecord = true;

            printf("[+] Ready to predict. Keep pressing keys and the program will guess which key was pressed\n");
            printf("    based on the captured audio from the microphone.\n");
            printf("[+] Predicting\n");
        }

        if (doRecord) {
            doRecord = false;
            audioLogger.record(kBufferSize_s);
        }
    };

    init();
    while (true) {
        update();
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

    return 0;
}
