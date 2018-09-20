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

// constants

constexpr float kBufferSize_s = 0.075f;
constexpr uint64_t kSampleRate = 96000;

constexpr uint64_t kRingBufferSize = 16*1024;
constexpr int bkgrStep_samples = 7;
constexpr int keyDuration_samples = 0.050f*kSampleRate;

constexpr uint64_t kBufferSize_frames = 2*AudioLogger::getBufferSize_frames(kSampleRate, kBufferSize_s) - 1;

constexpr auto kSamplesPerFrame = AudioLogger::kSamplesPerFrame;
constexpr auto kSamplesPerWaveform = kSamplesPerFrame*kBufferSize_frames;

// types

using TKey = int;
using TKeyWaveform = std::array<AudioLogger::Sample, kSamplesPerFrame*kBufferSize_frames>;
using TKeyHistory = std::vector<TKeyWaveform>;

// helpers

using TSum = double;
using TSum2 = double;
std::tuple<TSum, TSum2> calcSum(const TKeyWaveform & waveform, int is0, int is1) {
    TSum sum = 0.0f;
    TSum2 sum2 = 0.0f;
    for (int is = is0; is < is1; ++is) {
        auto a0 = waveform[is];
        sum += a0;
        sum2 += a0*a0;
    }

    return { sum, sum2 };
}

using TValueCC = double;
using TOffset = int;
TValueCC calcCC(
    const TKeyWaveform & waveform0,
    const TKeyWaveform & waveform1,
    TSum sum0, TSum2 sum02,
    int is00, int is0, int is1) {
    TValueCC cc = -1.0f;

    TSum sum1 = 0.0f;
    TSum2 sum12 = 0.0f;
    TSum2 sum01 = 0.0f;
    for (int is = 0; is < is1 - is0; ++is) {
        auto a0 = waveform0[is00 + is];
        auto a1 = waveform1[is0 + is];

        sum1 += a1;
        sum12 += a1*a1;
        sum01 += a0*a1;
    }

    int ncc = (is1 - is0);
    {
        double nom = sum01*ncc - sum0*sum1;
        double den2a = sum02*ncc - sum0*sum0;
        double den2b = sum12*ncc - sum1*sum1;
        cc = (nom)/(sqrt(den2a*den2b));
    }

    return cc;
}

std::tuple<TValueCC, TOffset> findBestCC(
    const TKeyWaveform & waveform0,
    const TKeyWaveform & waveform1,
    int is0, int is1,
    int alignWindow) {
    TOffset besto = -1;
    TValueCC bestcc = -1.0f;

    int is00 = waveform0.size()/2 - (is1 - is0)/2;
    auto [sum0, sum02] = calcSum(waveform0, is00, is00 + is1 - is0);

    int nWorkers = std::thread::hardware_concurrency();
    std::mutex mutex;
    std::vector<std::thread> workers(nWorkers);
    for (int i = 0; i < workers.size(); ++i) {
        auto & worker = workers[i];
        worker = std::thread([&, sum0 = sum0, sum02 = sum02, i]() {
            TOffset cbesto = -1;
            TValueCC cbestcc = -1.0f;

            for (int o = -alignWindow + i; o < alignWindow; o += nWorkers) {
                auto cc = calcCC(waveform0, waveform1, sum0, sum02, is00, is0 + o, is1 + o);
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

    return { bestcc, besto };
}

// globals

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

    TKey keyPressed = -1;
    std::map<TKey, TKeyHistory> keySoundHistoryAmpl;
    std::map<TKey, TKeyWaveform> keySoundAverageAmpl;

    bool doRecord = false;
    bool printStatus = true;
    bool isReadyToPredict = false;
    bool processingInput = true;

    int curFile = 0;

    // ring buffer
    int rbBegin = 0;
    float rbAverage = 0.0f;
    std::array<float, kRingBufferSize> rbSamples;
    rbSamples.fill(0.0f);

    AudioLogger audioLogger;

    struct WorkData {
        TKeyWaveform ampl;
        std::vector<int> positionsToPredict;
    };

    std::mutex mutex;
    std::deque<WorkData> workQueue;
    std::thread worker([&]() {
        int lastkey = -1;
        double lastcc = -1.0f;

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
                const auto & ampl = workData.ampl;
                const auto & positionsToPredict = workData.positionsToPredict;

                int nFramesPerWaveform = kBufferSize_frames;
                int alignWindow = 2*kSamplesPerFrame;

                for (int ipos = 0; ipos < positionsToPredict.size() ; ++ipos) {
                    int scmp0 = positionsToPredict[ipos] - kSamplesPerFrame/2;
                    int scmp1 = positionsToPredict[ipos] + kSamplesPerFrame/2;

                    char res = -1;
                    double maxcc = -1.0f;
                    for (const auto & ka : keySoundAverageAmpl) {
                        auto [bestcc, bestoffset] = findBestCC(keySoundAverageAmpl[ka.first], ampl, scmp0, scmp1, alignWindow);
                        if (bestcc > maxcc) {
                            res = ka.first;
                            maxcc = bestcc;
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
        if (frames.size() != kBufferSize_frames) {
            printf("Unexpected number of frames - %d, expected - %d. Should never happen\n",
                   (int) frames.size(), (int) kBufferSize_frames);
            return;
        }

        const int nFrames = frames.size();

        if (isReadyToPredict) {

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
                            positionsToPredict.push_back(f*kSamplesPerFrame + s);
                            //printf("Key press detected\n");
                        }
                    }
                }

                //printf("Average = %10.8f, max = %10.8f\n", rbAverage, amax);
            }

            if (positionsToPredict.size() > 0) {
                WorkData workData;
                auto & ampl = workData.ampl;
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
            history.push_back(TKeyWaveform());
            auto & ampl = history.back();
            for (int k = 0; k < nFrames; ++k) {
                std::copy(frames[k].begin(), frames[k].end(), ampl.begin() + k*kSamplesPerFrame);
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
                int nFramesPerWaveform = kBufferSize_frames;

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

                    for (int icur = 0; icur < kSamplesPerWaveform; ++icur) {
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
                        sum += p;
                        sum2 += p*p;
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

                int centerSample = kSamplesPerWaveform/2;

                printf("    - Centering waveforms at sample %d\n", centerSample);
                for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                    int offset = peakUsed[iwaveform] - centerSample;
                    //printf("        Offset for waveform %-4d = %-4d\n", iwaveform, offset);

                    auto newWaveform = TKeyWaveform();
                    auto & waveform = history[iwaveform];
                    for (int icur = 0; icur < kSamplesPerWaveform; ++icur) {
                        int iorg = icur + offset;

                        if (iorg >= 0 && iorg < kSamplesPerWaveform) {
                            newWaveform[icur] = waveform[iorg];
                        } else {
                            newWaveform[icur] = 0.0f;
                        }
                    }

                    waveform = std::move(newWaveform);
                }

                int alignWindow = centerSample/2;
                printf("    - Calculating CC pairs\n");
                printf("      Align window = %d\n", alignWindow);

                int bestw = -1;
                double bestccsum2 = -1.0f;
                std::map<int, std::map<int, std::tuple<TValueCC, TOffset>>> ccs;

                for (int alignToWaveform = 0; alignToWaveform < nWaveforms; ++alignToWaveform) {
                    ccs[alignToWaveform][alignToWaveform] = { 1.0f, 0 };

                    int is0 = centerSample - kSamplesPerFrame/2;
                    int is1 = centerSample + kSamplesPerFrame/2;

                    const auto & waveform0 = history[alignToWaveform];
                    auto [sum0, sum02] = calcSum(waveform0, is0, is1);

                    for (int iwaveform = alignToWaveform + 1; iwaveform < nWaveforms; ++iwaveform) {
                        const auto & waveform1 = history[iwaveform];
                        auto [bestcc, bestoffset] = findBestCC(waveform0, waveform1, is0, is1, alignWindow);

                        ccs[iwaveform][alignToWaveform] = { bestcc, bestoffset };
                        ccs[alignToWaveform][iwaveform] = { bestcc, -bestoffset };
                    }

                    double curccsum2 = 0.0f;
                    for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                        auto [cc, offset] = ccs[iwaveform][alignToWaveform];
                        curccsum2 += cc;
                    }

                    if (curccsum2 > bestccsum2) {
                        bestw = alignToWaveform;
                        bestccsum2 = curccsum2;
                    }
                }
                bestccsum2 /= nWaveforms;

                printf("    - Aligning all waveforms to waveform %d, (cost = %g)\n", bestw, bestccsum2);
                //std::ofstream fout(std::string("waveform_one_") + std::to_string(key) + ".plot");
                //for (auto & v : history[bestw]) fout << v << std::endl;
                //fout << std::endl;
                for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                    if (iwaveform == bestw) continue;

                    auto & waveform1 = history[iwaveform];
                    auto [cc, offset] = ccs[iwaveform][bestw];

                    auto newWaveform = TKeyWaveform();
                    for (int icur = 0; icur < kSamplesPerWaveform; ++icur) {
                        int iorg = icur + offset;

                        if (iorg >= 0 && iorg < kSamplesPerWaveform) {
                            newWaveform[icur] = waveform1[iorg];
                        } else {
                            newWaveform[icur] = 0.0f;
                        }
                    }

                    waveform1 = std::move(newWaveform);
                    //for (auto & v : waveform1) fout << v << std::endl;
                    //fout << std::endl;
                }

                printf("    - Calculating average waveform\n");
                double ccsum = 0.0f;
                double norm = 1.0f;
                auto & avgWaveform = keySoundAverageAmpl[key];
                avgWaveform.fill(0.0f);
                for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                    auto [cc, offset] = ccs[iwaveform][bestw];
                    if (cc < 0.4f) continue;
                    printf("        Adding waveform %d - cc = %g, offset = %d\n", iwaveform, cc, offset);
                    ccsum += cc*cc;
                    norm += cc*cc;
                    auto & waveform = history[iwaveform];
                    for (int is = 0; is < kSamplesPerWaveform; ++is) {
                        avgWaveform[is] += cc*cc*waveform[is];
                    }
                }

                norm = 1.0f/(norm);
                for (int is = 0; is < kSamplesPerWaveform; ++is) {
                    avgWaveform[is] *= norm;
                }

                if (ccsum*norm < 0.50f || (1.0f/norm < nWaveforms/3)) {
                    failedToTrain.push_back(key);
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
