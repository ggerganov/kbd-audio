/*! \file keytap.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "constants.h"
#include "common.h"
#include "audio_logger.h"

#include <map>
#include <mutex>
#include <cmath>
#include <string>
#include <chrono>
#include <thread>
#include <vector>
#include <deque>
#include <fstream>

//#define MY_DEBUG
//#define OUTPUT_WAVEFORMS

// globals

static int g_predictedKey = -1;
static bool g_isInitialized = false;

static std::function<int()> g_init;
static std::function<void()> g_update;
static std::function<bool()> g_mainUpdate;
static std::function<void(int)> g_handleKey;

int init() {
    if (g_isInitialized) return 1;

    return g_init();
}

void update() {
    if (g_isInitialized == false) return;

    g_update();
}

void mainUpdate() {
    g_mainUpdate();
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
    printf("Usage: %s input.kbd [input2.kbd ...] [-cN] [-pF] [-tF]\n", argv[0]);
    printf("    -cN - select capture device N\n");
    printf("    -CN - select number N of capture channels to use\n");
    printf("    -pF - prediction threshold: CC > F\n");
    printf("    -tF - background threshold: ampl > F*avg_background\n");
    printf("\n");

    if (argc < 2) {
        return -127;
    }

    auto argm = parseCmdArguments(argc, argv);
    int captureId = argm["c"].empty() ? 0 : std::stoi(argm["c"]);
    int nChannels = argm["C"].empty() ? 0 : std::stoi(argm["C"]);

    std::map<int, std::ifstream> fins;
    for (int i = 0; i < argc - 1; ++i) {
        if (argv[i + 1][0] == '-') continue;
        printf("Opening file '%s'\n", argv[i + 1]);
        fins[i] = std::ifstream(argv[i + 1], std::ios::binary);
        if (fins[i].good() == false) {
            printf("Failed to open input file: '%s'\n", argv[i + 1]);
            return -2;
        }

        {
            int bufferSize_frames = 1;
            fins[i].read((char *)(&bufferSize_frames), sizeof(bufferSize_frames));
            if (bufferSize_frames != kBufferSizeTrain_frames) {
                printf("Buffer size in file (%d) does not match the expected one (%d)\n", bufferSize_frames, (int) kBufferSizeTrain_frames);
                return -1;
            }
        }
    }

    TKey keyPressed = -1;
    std::map<TKey, TKeyHistoryF> keySoundHistoryAmpl;
    std::map<TKey, TKeyWaveformF> keySoundAverageAmpl;

    int ntest = 0;

    bool doRecord = false;
    bool isReadyToPredict = false;
    bool finishApp = false;
    bool processingInput = true;

    int curFile = 0;

    float amplMin = 0.0f;
    float amplMax = 0.0f;
    float thresholdCC = argm["p"].empty() ? 0.35f : std::stof(argm["p"]);
    float thresholdBackground = argm["t"].empty() ? 10.0f : std::stof(argm["t"]);

    // ring buffer
    int rbBegin = 0;
    float rbAverage = 0.0f;
    std::array<float, kBkgrRingBufferSize> rbSamples;
    rbSamples.fill(0.0f);

    // Train data
    bool isAcquiringTrainData = false;
    std::map<int, int> nTimes;
    size_t totalSize_bytes = 0;
    std::ofstream foutTrain("train_default.kbd", std::ios::binary);
    {
        int x = kBufferSizeTrain_frames;
        foutTrain.write((char *)(&x), sizeof(x));
    }

    AudioLogger audioLogger;

    struct WorkData {
        TKeyWaveformF ampl;
        std::vector<int> positionsToPredict;
    };

    std::mutex mutex;
    std::deque<WorkData> workQueue;
    std::thread worker([&]() {
        int lastkey = -1;
        double lastcc = -1.0f;

        while (finishApp == false) {
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
                const auto & ampl = workData.ampl;
                const auto & positionsToPredict = workData.positionsToPredict;

                //int alignWindow = kSamplesPerFrame/2;
                int alignWindow = 64;

                for (int ipos = 0; ipos < (int) positionsToPredict.size(); ++ipos) {
                    auto curPos = positionsToPredict[ipos];
                    int scmp0 = curPos - kSamplesPerFrame;
                    int scmp1 = curPos + kSamplesPerFrame;

                    char res = -1;
                    TValueCC maxcc = -1.0f;
                    //TOffset offs = 0;
                    TKeyConfidenceMap keyConfidenceTmp;
                    for (const auto & ka : keySoundAverageAmpl) {
                        //auto [bestcc, bestoffset] = findBestCC(keySoundAverageAmpl[ka.first], ampl, scmp0, scmp1, alignWindow);
                        auto ret = findBestCC(keySoundAverageAmpl[ka.first], ampl, scmp0, scmp1, alignWindow);
                        auto bestcc     = std::get<0>(ret);
                        //auto bestoffset = std::get<1>(ret);

                        if (bestcc > maxcc) {
                            res = ka.first;
                            maxcc = bestcc;
                            //offs = bestoffset;
                        }
                        keyConfidenceTmp[ka.first] = bestcc;
                    }

                    if (maxcc > thresholdCC) {
                        if (lastkey != res || lastcc != maxcc) {
                            printf("    Prediction: '%c'        (%8.5g), ntest = %d\n", res, maxcc, ntest);
                        }
                        lastkey = res;
                        lastcc = maxcc;
                    }
                    ++ntest;
                }

            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
    });

    AudioLogger::Callback cbAudio = [&](const AudioLogger::Record & frames) {
        if (isAcquiringTrainData) {
            foutTrain.write((char *)(&keyPressed), sizeof(keyPressed));
            for (const auto & frame : frames) {
                totalSize_bytes += sizeof(frame[0])*frame.size();
                foutTrain.write((char *)(frame.data()), sizeof(frame[0])*frame.size());
                foutTrain.flush();
            }
            ++nTimes[keyPressed];

            printf("Last recorded key - %3d '%s'. Total times recorded so far - %3d. Total data saved: %g MB\n",
                   keyPressed, kKeyText.at(keyPressed), nTimes[keyPressed], ((float)(totalSize_bytes)/1024.0f/1024.0f));

            keyPressed = -1;
            return;
        }

        if (frames.size() != kBufferSizeTrain_frames && isReadyToPredict == false) {
            printf("Unexpected number of frames - %d, expected - %d. Should never happen\n",
                   (int) frames.size(), (int) kBufferSizeTrain_frames);
            return;
        }

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

                int nFrames = frames.size();

                auto _acc = [](const AudioLogger::Record & r, int id) { return std::abs(r[id/kSamplesPerFrame][id%kSamplesPerFrame]); };

                int k = kSamplesPerFrame;
                std::deque<int> que(k);
                for (int i = 0; i < nFrames*kSamplesPerFrame; ++i) {
                    if (i < k) {
                        while((!que.empty()) && _acc(frames, i) >= _acc(frames, que.back())) {
                            que.pop_back();
                        }
                        que.push_back(i);
                    } else {
                        while((!que.empty()) && que.front() <= i - k) {
                            que.pop_front();
                        }

                        while((!que.empty()) && _acc(frames, i) >= _acc(frames, que.back())) {
                            que.pop_back();
                        }

                        que.push_back(i);

                        int itest = i - k/2;
                        if (itest >= (0.5*kSamplesPerWaveformTrain - kSamplesPerFrame) && itest < (0.5*kSamplesPerWaveformTrain + kSamplesPerFrame) && que.front() == itest) {
                            auto acur = _acc(frames, itest);
                            if (acur > thresholdBackground*rbAverage){
                                positionsToPredict.push_back(itest);
                            }
                        }
                    }
                }
            }

            if (positionsToPredict.size() > 0) {
                WorkData workData;
                auto & ampl = workData.ampl;
                ampl.resize(nFrames*kSamplesPerFrame);
                for (int k = 0; k < nFrames; ++k) {
                    std::copy(frames[k].begin(), frames[k].end(), ampl.begin() + k*kSamplesPerFrame);
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
            history.push_back(TKeyWaveformF());
            auto & ampl = history.back();
            ampl.resize(nFrames*kSamplesPerFrame);
            for (int k = 0; k < nFrames; ++k) {
                std::copy(frames[k].begin(), frames[k].end(), ampl.begin() + k*kSamplesPerFrame);
            }
        }

        keyPressed = -1;
    };

    g_init = [&]() {
        AudioLogger::Parameters parameters;
        parameters.callback = cbAudio;
        parameters.captureId = captureId;
        parameters.nChannels = nChannels;
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

    g_update = [&]() {
        if (isAcquiringTrainData) {
            return;
        }

        if (processingInput) {
            if (keyPressed == -1) {
                AudioLogger::Frame frame;
                AudioLogger::Record record;
                fins[curFile].read((char *)(&keyPressed), sizeof(keyPressed));
                if (fins[curFile].eof()) {
                    ++curFile;
                    if (curFile >= (int) fins.size()) {
                        processingInput = false;
                    }
                } else {
                    printf("%c", keyPressed);
                    fflush(stdout);
                    for (int i = 0; i < kBufferSizeTrain_frames; ++i) {
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
                int nFramesPerWaveform = kBufferSizeTrain_frames;

                printf("    - Training key '%c'\n", key);
                printf("    - History size = %d key waveforms\n", nWaveforms);
                printf("    - Frames per key waveform   = %d\n", nFramesPerWaveform);
                printf("    - Total frames available    = %d\n", nWaveforms*nFramesPerWaveform);
                printf("    - Samples per frame         = %d\n", (int) kSamplesPerFrame);
                printf("    - Total samples available   = %d\n", (int) (nWaveforms*nFramesPerWaveform*kSamplesPerFrame));

                printf("    - Estimating waveform peaks ...\n");
                std::vector<int> peakSum;
                std::vector<int> peakMax;

                peakSum.clear();
                peakMax.clear();
                for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                    int isum = -1;
                    double asum = 0.0f;
                    double aisum = 0.0f;

                    int imax = -1;
                    double amax = 0.0f;

                    const auto & waveform = history[iwaveform];

                    for (int icur = 0; icur < kSamplesPerWaveformTrain; ++icur) {
                        double acur = std::abs(waveform[icur]);
                        double acur2 = acur*acur;

                        asum += acur2;
                        aisum += acur2*icur;

                        if (acur > amax) {
                            amax = acur;
                            imax = icur;
                        }
                    }

                    isum = aisum/asum;

                    peakSum.push_back(isum);
                    peakMax.push_back(imax);
                    //printf("        Estimated peak: %d (method - sum), %d (method - max)\n", isum, imax);
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

                double stdevSum = calcStdev(peakSum);
                double stdevMax = calcStdev(peakMax);

                printf("    - Stdev of estimated peaks: %g (sum) vs %g (max)\n", stdevSum, stdevMax);

                const auto & peakUsed = peakMax;
                printf("    - Using 'max' estimation\n");

                int centerSample = kSamplesPerWaveformTrain/2;

                printf("    - Centering waveforms at sample %d\n", centerSample);
                for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                    int offset = peakUsed[iwaveform] - centerSample;
                    //printf("        Offset for waveform %-4d = %-4d\n", iwaveform, offset);

                    auto newWaveform = TKeyWaveformF();
                    newWaveform.resize(kSamplesPerWaveformTrain);
                    auto & waveform = history[iwaveform];
                    for (int icur = 0; icur < kSamplesPerWaveformTrain; ++icur) {
                        int iorg = icur + offset;

                        if (iorg >= 0 && iorg < kSamplesPerWaveformTrain) {
                            newWaveform[icur] = waveform[iorg];
                        } else {
                            newWaveform[icur] = 0.0f;
                        }
                    }

                    waveform = std::move(newWaveform);
                }

                int alignWindow = 64;
                printf("    - Calculating CC pairs\n");
                printf("      Align window = %d\n", alignWindow);

                int bestw = -1;
                int ntrain = 0;
                double bestccsum = -1.0f;
                //double bestosum = 1e10;
                std::map<int, std::map<int, std::tuple<TValueCC, TOffset>>> ccs;

                for (int alignToWaveform = 0; alignToWaveform < nWaveforms; ++alignToWaveform) {
                    ccs[alignToWaveform][alignToWaveform] = std::tuple<TValueCC, TOffset>(1.0f, 0);

                    int is0 = centerSample - kSamplesPerFrame;
                    int is1 = centerSample + kSamplesPerFrame;

                    const auto & waveform0 = history[alignToWaveform];

                    for (int iwaveform = alignToWaveform + 1; iwaveform < nWaveforms; ++iwaveform) {
                        const auto & waveform1 = history[iwaveform];
                        //auto [bestcc, bestoffset] = findBestCC(waveform0, waveform1, is0, is1, alignWindow);
                        auto ret = findBestCC(waveform0, waveform1, is0, is1, alignWindow);
                        auto bestcc     = std::get<0>(ret);
                        auto bestoffset = std::get<1>(ret);

                        ccs[iwaveform][alignToWaveform] = std::tuple<TValueCC, TOffset>(bestcc, bestoffset);
                        ccs[alignToWaveform][iwaveform] = std::tuple<TValueCC, TOffset>(bestcc, -bestoffset);
                    }

                    int curntrain = 0;
                    double curccsum = 0.0;
                    double curosum = 0.0;
                    for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                        //auto [cc, offset] = ccs[iwaveform][alignToWaveform];
                        auto cc     = std::get<0>(ccs[iwaveform][alignToWaveform]);
                        auto offset = std::get<1>(ccs[iwaveform][alignToWaveform]);

                        if (std::abs(offset) > 50) continue;
                        ++curntrain;
                        curccsum += cc*cc;
                        curosum += offset*offset;
                    }

                    if (curccsum > bestccsum) {
                    //if (curosum < bestosum) {
                        ntrain = curntrain;
                        bestw = alignToWaveform;
                        bestccsum = curccsum;
                        //bestosum = curosum;
                    }
                }
                bestccsum = sqrt(bestccsum/ntrain);

                printf("    - Aligning all waveforms to waveform %d, (cost = %g)\n", bestw, bestccsum);
#ifdef OUTPUT_WAVEFORMS
                std::ofstream fout(std::string("waveform_one_") + std::to_string(key) + ".plot");
                for (auto & v : history[bestw]) fout << v << std::endl;
                fout << std::endl;
#endif
                for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                    if (iwaveform == bestw) continue;

                    auto & waveform1 = history[iwaveform];
                    //auto [cc, offset] = ccs[iwaveform][bestw];
                    //auto cc     = std::get<0>(ccs[iwaveform][bestw]);
                    auto offset = std::get<1>(ccs[iwaveform][bestw]);

                    auto newWaveform = TKeyWaveformF();
                    newWaveform.resize(kSamplesPerWaveformTrain);
                    for (int icur = 0; icur < kSamplesPerWaveformTrain; ++icur) {
                        int iorg = icur + offset;

                        if (iorg >= 0 && iorg < kSamplesPerWaveformTrain) {
                            newWaveform[icur] = waveform1[iorg];
                        } else {
                            newWaveform[icur] = 0.0f;
                        }
                    }

                    waveform1 = std::move(newWaveform);
#ifdef OUTPUT_WAVEFORMS
                    for (auto & v : waveform1) fout << v << std::endl;
                    fout << std::endl;
#endif
                }

                printf("    - Calculating average waveform\n");
                double ccsum = 0.0f;
                double norm = 0.0f;
                auto & avgWaveform = keySoundAverageAmpl[key];
                avgWaveform.resize(kSamplesPerWaveformTrain);
                std::fill(avgWaveform.begin(), avgWaveform.end(), 0.0f);
                for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                    //auto [cc, offset] = ccs[iwaveform][bestw];
                    auto cc     = std::get<0>(ccs[iwaveform][bestw]);
                    auto offset = std::get<1>(ccs[iwaveform][bestw]);

                    //if (std::abs(offset) > 5) continue;
                    printf("        Adding waveform %d - cc = %g, offset = %ld\n", iwaveform, cc, offset);
                    ccsum += cc*cc;
                    norm += cc*cc;
                    auto & waveform = history[iwaveform];
                    for (int is = 0; is < kSamplesPerWaveformTrain; ++is) {
                        avgWaveform[is] += cc*cc*waveform[is];
                    }
                }

                norm = 1.0f/(norm);
                for (int is = 0; is < kSamplesPerWaveformTrain; ++is) {
                    avgWaveform[is] *= norm;
                    if (avgWaveform[is] > amplMax) amplMax = avgWaveform[is];
                    if (avgWaveform[is] < amplMin) amplMin = avgWaveform[is];
                }

#ifdef OUTPUT_WAVEFORMS
                {
                    std::ofstream fout(std::string("waveform_avg_") + std::to_string(key) + ".plot");
                    for (auto & v : avgWaveform) fout << v << std::endl;
                }
#endif

                if (ccsum*norm < 0.50f || (1.0f/norm < nWaveforms/3)) {
                    failedToTrain.push_back(key);
                }

                printf("\n");
            };

            for (const auto & kh : keySoundHistoryAmpl) {
                if (kh.second.size() > 2) {
                    trainKey(kh.first);
                } else {
                    printf("[!] Key '%s' does not have enough training data. Need at least 3 presses\n", kKeyText.at(kh.first));
                    failedToTrain.push_back(kh.first);
                }
            }
            printf("Failed to train the following keys: ");
            for (auto & k : failedToTrain) printf("'%c' ", k);
            printf("\n");
            isReadyToPredict = true;
            doRecord = true;

            amplMax = std::max(amplMax, -amplMin);
            amplMin = -std::max(amplMax, -amplMin);

            for (auto & kh : keySoundAverageAmpl) {
                float curAmplMax = 0.0f;
                for (const auto & v : kh.second) if (std::abs(v) > curAmplMax) curAmplMax = std::abs(v);
                for (auto & v : kh.second) v = (v/curAmplMax)*amplMax;
            }

            audioLogger.resume();

            printf("[+] Ready to predict. Keep pressing keys and the program will guess which key was pressed\n");
            printf("    based on the captured audio from the microphone.\n");
            printf("[+] Predicting\n");
        }

        if (doRecord) {
            doRecord = false;
            audioLogger.record(kBufferSizeTrain_s, getBufferSize_frames(kSampleRate, kBufferSizeTrain_s) - 1);
        }
    };

    g_mainUpdate = [&]() {
        if (finishApp) return false;

        update();

        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        return true;
    };

    init();
#ifdef __EMSCRIPTEN__
    emscripten_set_main_loop(mainUpdate, 60, 1);
#else
    while (true) {
        if (g_mainUpdate() == false) break;
    }
#endif

    worker.join();

    printf("[+] Terminated");

    return 0;
}
