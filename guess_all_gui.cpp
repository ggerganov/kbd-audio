/*! \file guess_all_gui.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "constants.h"
#include "audio_logger.h"

#include "imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl3.h"

#include <SDL.h>
#include <GL/gl3w.h>

#include <map>
#include <cmath>
#include <string>
#include <chrono>
#include <thread>
#include <vector>
#include <deque>
#include <fstream>

#define MY_DEBUG
#define OUTPUT_WAVEFORMS

// constants

constexpr float kBufferSize_s = 0.075f;
constexpr uint64_t kSampleRate = 24000;

constexpr uint64_t kRingBufferSize = 128*1024;
constexpr int bkgrStep_samples = 7;
constexpr int keyDuration_samples = 0.005f*kSampleRate;

constexpr uint64_t kBufferSize_frames = 2*AudioLogger::getBufferSize_frames(kSampleRate, kBufferSize_s) - 1;

constexpr auto kSamplesPerFrame = AudioLogger::kSamplesPerFrame;
constexpr auto kSamplesPerWaveform = kSamplesPerFrame*kBufferSize_frames;

static const std::vector<float> kRowOffset = { 0.0f, 1.5f, 1.8f, 2.1f, 5.5f };
static const std::vector<std::vector<int>> kKeyboard = {
    { '`', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '-', '=', 127},
    { 'q', 'w', 'e', 'r', 't', 'y', 'u', 'i', 'o', 'p', '[', ']', '\\', },
    { 'a', 's', 'd', 'f', 'g', 'h', 'j', 'k', 'l', ';', '\'', '\n' },
    { 'z', 'x', 'c', 'v', 'b', 'n', 'm', ',', '.', '/', },
    { ' ', },
};

// types

using TConfidence = float;

using TKey = int;
using TKeyWaveform = std::vector<AudioLogger::Sample>;
using TKeyHistory = std::vector<TKeyWaveform>;
using TKeyConfidenceMap = std::map<TKey, TConfidence>;

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

#ifdef MY_DEBUG
        if (is00 + is < 0 || is00 + is >= waveform0.size()) printf("BUG 0\n");
        if (is0 + is < 0 || is0 + is >= waveform1.size()) {
            printf("BUG 1\n");
            printf("%d %d %d\n", is0, is, (int) waveform1.size());
        }
#endif

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
}

int main(int argc, char ** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.kbd [input2.kbd ...]\n", argv[0]);
        return -127;
    }

    if (SDL_Init(SDL_INIT_VIDEO|SDL_INIT_TIMER) != 0) {
        printf("Error: %s\n", SDL_GetError());
        return -1;
    }

    int windowSizeX = 1280;
    int windowSizeY = 800;

#if __APPLE__
    // GL 3.2 Core + GLSL 150
    const char* glsl_version = "#version 150";
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);
#else
    // GL 3.0 + GLSL 130
    const char* glsl_version = "#version 130";
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, 0);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
#endif

    // Create window with graphics context
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);
    SDL_DisplayMode current;
    SDL_GetCurrentDisplayMode(0, &current);
    SDL_Window* window = SDL_CreateWindow("Keytap", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, windowSizeX, windowSizeY, SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE);
    SDL_GLContext gl_context = SDL_GL_CreateContext(window);
    SDL_GL_SetSwapInterval(1); // Enable vsync

    // Initialize OpenGL loader
    bool err = gl3wInit() != 0;
    if (err) {
        fprintf(stderr, "Failed to initialize OpenGL loader!\n");
        return 1;
    }

    // Setup Dear ImGui binding
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard Controls

    ImGui_ImplSDL2_InitForOpenGL(window, gl_context);
    ImGui_ImplOpenGL3_Init(glsl_version);

    // Setup style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsClassic();

    ImFontConfig fontConfig;
    fontConfig.SizePixels = 24.0f;
    ImGui::GetIO().Fonts->AddFontDefault(&fontConfig);

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
    TKeyConfidenceMap keyConfidence;
    std::map<TKey, TKeyHistory> keySoundHistoryAmpl;
    std::map<TKey, TKeyWaveform> keySoundAverageAmpl;

    bool doRecord = false;
    bool printStatus = true;
    bool isReadyToPredict = false;
    bool processingInput = true;

    int curFile = 0;

    int predictedKey = -1;
    float thresholdCC = 0.5f;
    float thresholdBackground = 10.0f;
    TValueCC predictedCC = -1.0f;
    auto tLastDetectedKeyStroke = std::chrono::high_resolution_clock::now();

    // ring buffer
    int rbBegin = 0;
    double rbAverage = 0.0f;
    std::array<double, kRingBufferSize> rbSamples;
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
                int alignWindow = ((nFramesPerWaveform/2)*kSamplesPerFrame)/2;

                for (int ipos = 0; ipos < positionsToPredict.size() ; ++ipos) {
                    int scmp0 = positionsToPredict[ipos] - kSamplesPerFrame/2;
                    int scmp1 = positionsToPredict[ipos] + kSamplesPerFrame/2;

                    char res = -1;
                    double maxcc = -1.0f;
                    TKeyConfidenceMap keyConfidenceTmp;
                    for (const auto & ka : keySoundAverageAmpl) {
                        auto [bestcc, bestoffset] = findBestCC(keySoundAverageAmpl[ka.first], ampl, scmp0, scmp1, alignWindow);
                        //printf(" %8.4f ", bestcc);
                        if (bestcc > maxcc) {
                            res = ka.first;
                            maxcc = bestcc;
                        }
                        keyConfidenceTmp[ka.first] = bestcc;
                    }
                    //printf("\n");

                    if (maxcc > thresholdCC) {
                        if (lastkey != res || lastcc != maxcc) {
                            printf("    Prediction: '%c'        (%8.5g)\n", res, maxcc);
                            predictedKey = res;
                            predictedCC = maxcc;
                            for (auto & c : keyConfidenceTmp) {
                                keyConfidence[c.first] = std::pow(c.second/maxcc, 4.0);
                            }
                        }
                        lastkey = res;
                        lastcc = maxcc;
                    }
                }
            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
    });

    AudioLogger::Callback cbAudio = [&](const AudioLogger::Record & frames) {
        if (frames.size() != kBufferSize_frames && isReadyToPredict == false) {
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
                //for (int f = nFrames2 - nFrames2/2; f <= nFrames2 + nFrames2/2; ++f) {
                for (int f = 2; f < frames.size() - 2; ++f) {
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
                        if (acur > thresholdBackground*rbAverage && acur > 0.5f*amax &&
                            ((s == frames[f].size() - 1) || (acur > frames[f][s-1] && acur > frames[f][s+1]))) {
                            skip_samples = keyDuration_samples;
                            positionsToPredict.push_back(f*kSamplesPerFrame + s);
                            tLastDetectedKeyStroke = std::chrono::high_resolution_clock::now();
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
            history.push_back(TKeyWaveform());
            auto & ampl = history.back();
            ampl.resize(nFrames*kSamplesPerFrame);
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
            predictedKey = -1;
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
                    newWaveform.resize(kSamplesPerWaveform);
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

                    for (int iwaveform = alignToWaveform + 1; iwaveform < nWaveforms; ++iwaveform) {
                        const auto & waveform1 = history[iwaveform];
                        auto [bestcc, bestoffset] = findBestCC(waveform0, waveform1, is0, is1, alignWindow);

                        ccs[iwaveform][alignToWaveform] = { bestcc, bestoffset };
                        ccs[alignToWaveform][iwaveform] = { bestcc, -bestoffset };
                    }

                    double curccsum2 = 0.0f;
                    for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                        auto [cc, offset] = ccs[iwaveform][alignToWaveform];
                        if (std::abs(offset) > 5) continue;
                        curccsum2 += cc;
                    }

                    if (curccsum2 > bestccsum2) {
                        bestw = alignToWaveform;
                        bestccsum2 = curccsum2;
                    }
                }
                bestccsum2 /= nWaveforms;

                printf("    - Aligning all waveforms to waveform %d, (cost = %g)\n", bestw, bestccsum2);
#ifdef OUTPUT_WAVEFORMS
                std::ofstream fout(std::string("waveform_one_") + std::to_string(key) + ".plot");
                for (auto & v : history[bestw]) fout << v << std::endl;
                fout << std::endl;
#endif
                for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                    if (iwaveform == bestw) continue;

                    auto & waveform1 = history[iwaveform];
                    auto [cc, offset] = ccs[iwaveform][bestw];

                    auto newWaveform = TKeyWaveform();
                    newWaveform.resize(kSamplesPerWaveform);
                    for (int icur = 0; icur < kSamplesPerWaveform; ++icur) {
                        int iorg = icur + offset;

                        if (iorg >= 0 && iorg < kSamplesPerWaveform) {
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
                double norm = 1.0f;
                auto & avgWaveform = keySoundAverageAmpl[key];
                avgWaveform.resize(kSamplesPerWaveform);
                std::fill(avgWaveform.begin(), avgWaveform.end(), 0.0f);
                for (int iwaveform = 0; iwaveform < nWaveforms; ++iwaveform) {
                    auto [cc, offset] = ccs[iwaveform][bestw];
                    if (cc < 0.50f || std::abs(offset) > 5) continue;
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
            audioLogger.recordSym(0.20f);
        }
    };

    bool done = false;

    init();
    while (done == false) {
        update();

        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            ImGui_ImplSDL2_ProcessEvent(&event);
            switch (event.type) {
                case SDL_KEYDOWN:
                    if (keyPressed == -1) {
                        //auto t1 = std::chrono::high_resolution_clock::now();
                        //printf("Event: %d\n", (int) std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());
                        keyPressed = event.key.keysym.sym;
                        audioLogger.record(kBufferSize_s);
                    }
                    break;
                case SDL_QUIT:
                    done = true;
                    break;
                case SDL_WINDOWEVENT:
                    if (event.window.event == SDL_WINDOWEVENT_CLOSE && event.window.windowID == SDL_GetWindowID(window)) done = true;
                    break;
            };
        }

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplSDL2_NewFrame(window);
        ImGui::NewFrame();

        ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Once);
        ImGui::SetNextWindowSize(ImVec2(windowSizeX, windowSizeY), ImGuiCond_Once);
        ImGui::Begin("Main");
        if (isReadyToPredict == false) {
            ImGui::Text("Training ... Please wait");
        } else {

            ImGui::Text("Last predicted key:       %s (%8.6f)\n", kKeyText.at(predictedKey), predictedCC);
            ImGui::SliderFloat("Threshold CC", &thresholdCC, 0.1f, 1.0f);
            auto tNow = std::chrono::high_resolution_clock::now();
            ImGui::Text("Last detected key stroke: %5.3f seconds ago\n",
                        (float)(std::chrono::duration_cast<std::chrono::milliseconds>(tNow - tLastDetectedKeyStroke).count()/1000.0f));
            ImGui::Text("Average background level: %16.13f\n", rbAverage);
            ImGui::SliderFloat("Threshold background", &thresholdBackground, 0.1f, 30.0f);
            ImGui::Text("\n");

            auto drawList = ImGui::GetWindowDrawList();

            auto p0 = ImGui::GetCursorScreenPos();
            auto p1 = p0;
            float ox = 0.0f;
            float oy = p0.y;
            float bx = 64.0f;
            float by = 64.0f;

            for (int rid = 0; rid < kKeyboard.size(); ++rid) {
                const auto & row = kKeyboard[rid];
                ox = p0.x + kRowOffset[rid]*bx;
                for (const auto & button : row) {
                    int key = button;
                    auto & conf = keyConfidence[key];
                    p1.x = ox;
                    p1.y = oy;
                    auto p2 = p1;
                    auto tw = ImGui::CalcTextSize(kKeyText.at(key));
                    p2.x += bx + tw.x;
                    p2.y += by;
                    ox += bx + tw.x;
                    drawList->AddRectFilled(p1, p2, ImGui::ColorConvertFloat4ToU32({ 0.2f, 1.0f, 0.2f, conf }), 0.0f, 0);
                    drawList->AddRect(p1, p2, ImGui::ColorConvertFloat4ToU32({ 1.0f, 1.0f, 1.0f, 1.0f }), 0.0f, 0, 2.0f);
                    ImGui::SetCursorScreenPos({0.5f*(p1.x + p2.x - tw.x), 0.5f*(p1.y + p2.y - tw.y)});
                    if (key == predictedKey) {
                        ImGui::TextColored({1.0f, 0.0f, 0.0f, 1.0f}, "%s", kKeyText.at(key));
                    } else {
                        ImGui::Text("%s", kKeyText.at(key));
                    }
                    conf *= 0.99f;
                }
                oy += by;
            }
        }

        ImGui::End();

        ImGui::Render();
        SDL_GL_MakeCurrent(window, gl_context);
        glViewport(0, 0, (int) io.DisplaySize.x, (int) io.DisplaySize.y);
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        SDL_GL_SwapWindow(window);
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplSDL2_Shutdown();
    ImGui::DestroyContext();

    SDL_GL_DeleteContext(gl_context);
    SDL_DestroyWindow(window);
    SDL_Quit();

    printf("Kill with Ctrl + C\n");
    while(true) { }

    return 0;
}
