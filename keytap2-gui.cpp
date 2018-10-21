/*! \file keytap2-gui.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl3.h"

#include <SDL.h>
#include <GL/gl3w.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <deque>
#include <fstream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#define MY_DEBUG

constexpr int kKeyPressWidth_samples = 256;

struct stMatch;
struct stKeyPressData;
struct stWaveformView;

using TSum                  = int64_t;
using TSum2                 = int64_t;
using TCC                   = double;
using TOffset               = int64_t;
using TMatch                = stMatch;
using TSimilarityMap        = std::vector<std::vector<TMatch>>;

using TClusterId            = int32_t;
using TSampleInput          = float;
using TSample               = int32_t;
using TWaveform             = std::vector<TSample>;
using TWaveformView         = stWaveformView;
using TKeyPressPosition     = int64_t;
using TKeyPressData         = stKeyPressData;
using TKeyPressCollection   = std::vector<TKeyPressData>;

struct stMatch {
    TCC     cc;
    TOffset offset;
};

struct stWaveformView {
    const TSample * samples;
    int64_t         n;
};

struct stKeyPressData {
    TWaveformView       waveform;
    TKeyPressPosition   pos;
    TCC                 ccAvg;
};

template <typename T>
float toSeconds(T t0, T t1) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()/1024.0f;
}

inline float frand() { return ((float)rand())/RAND_MAX; }

TWaveformView getView(const TWaveform & waveform, int64_t idx) { return { waveform.data() + idx, (int64_t) waveform.size() - idx }; }

TWaveformView getView(const TWaveform & waveform, int64_t idx, int64_t len) { return { waveform.data() + idx, len }; }

bool readFromFile(const std::string & fname, TWaveform & res) {
    std::ifstream fin(fname, std::ios::binary | std::ios::ate);
    if (fin.good() == false) {
        return false;
    }

    {
        std::streamsize size = fin.tellg();
        fin.seekg(0, std::ios::beg);

        static_assert(std::is_same<TSampleInput, float>::value, "TSampleInput not recognised");
        static_assert(
            std::is_same<TSample, float>::value
            || std::is_same<TSample, int16_t>::value
            || std::is_same<TSample, int32_t>::value
                      , "TSampleInput not recognised");

        if (std::is_same<TSample, int16_t>::value) {
            std::vector<TSampleInput> buf(size/sizeof(TSampleInput));
            res.resize(size/sizeof(TSampleInput));
            fin.read((char *)(buf.data()), size);
            double amax = 0.0f;
            for (auto i = 0; i < buf.size(); ++i) if (std::abs(buf[i]) > amax) amax = std::abs(buf[i]);
            for (auto i = 0; i < buf.size(); ++i) res[i] = std::round(32000.0*(buf[i]/amax));
            //double asum = 0.0f;
            //for (auto i = 0; i < buf.size(); ++i) asum += std::abs(buf[i]); asum /= buf.size(); asum *= 10.0;
            //for (auto i = 0; i < buf.size(); ++i) res[i] = std::round(2000.0*(buf[i]/asum));
        } else if (std::is_same<TSample, int32_t>::value) {
            std::vector<TSampleInput> buf(size/sizeof(TSampleInput));
            res.resize(size/sizeof(TSampleInput));
            fin.read((char *)(buf.data()), size);
            double amax = 0.0f;
            for (auto i = 0; i < buf.size(); ++i) if (std::abs(buf[i]) > amax) amax = std::abs(buf[i]);
            for (auto i = 0; i < buf.size(); ++i) res[i] = std::round(32000.0*(buf[i]/amax));
        } else if (std::is_same<TSample, float>::value) {
            res.resize(size/sizeof(TSample));
            fin.read((char *)(res.data()), size);
        } else {
        }
    }

    fin.close();

    return true;
}

bool findKeyPresses(const TWaveformView & waveform, TKeyPressCollection & res) {
    res.clear();

    int rbBegin = 0;
    double rbAverage = 0.0;
    std::array<double, 4*1024> rbSamples;
    rbSamples.fill(0.0);

    int k = 1024;
    double thresholdBackground = 10.0;

    std::deque<int64_t> que(k);
    auto [samples, n] = waveform;
    for (int64_t i = 0; i < n; ++i) {
        {
            rbAverage *= rbSamples.size();
            rbAverage -= rbSamples[rbBegin];
            double acur = std::abs(samples[i]);
            rbSamples[rbBegin] = acur;
            rbAverage += acur;
            rbAverage /= rbSamples.size();
            if (++rbBegin >= rbSamples.size()) {
                rbBegin = 0;
            }
        }

        if (i < k) {
            while((!que.empty()) && samples[i] >= samples[que.back()]) {
                que.pop_back();
            }
            que.push_back(i);
        } else {
            while((!que.empty()) && que.front() <= i - k) {
                que.pop_front();
            }

            while((!que.empty()) && samples[i] >= samples[que.back()]) {
                que.pop_back();
            }

            que.push_back(i);

            int64_t itest = i - k/2;
            if (itest >= 2*k && itest < n - 2*k && que.front() == itest) {
                double acur = samples[itest];
                if (acur > thresholdBackground*rbAverage){
                    TKeyPressData entry;
                    entry.waveform = waveform;
                    entry.pos = itest;
                    entry.ccAvg = 0.0;
                    res.emplace_back(std::move(entry));
                }
            }
        }
    }

    return true;
}

bool findKeyPresses(const TWaveform & waveform, TKeyPressCollection & res) {
    return findKeyPresses(getView(waveform, 0), res);
}

bool dumpKeyPresses(const std::string & fname, const TKeyPressCollection & data) {
    std::ofstream fout(fname);
    for (auto & k : data) {
        fout << k.pos << " 1" << std::endl;
    }
    fout.close();
    return true;
}

std::tuple<TSum, TSum2> calcSum(const TWaveformView & waveform) {
    TSum sum = 0.0f;
    TSum2 sum2 = 0.0f;
    auto [samples, n] = waveform;
    for (int64_t is = 0; is < n; ++is) {
        auto a0 = samples[is];
        sum += a0;
        sum2 += a0*a0;
    }

    return { sum, sum2 };
}

TCC calcCC(
    const TWaveformView & waveform0,
    const TWaveformView & waveform1,
    TSum sum0, TSum2 sum02) {
    TCC cc = -1.0f;

    TSum sum1 = 0.0f;
    TSum2 sum12 = 0.0f;
    TSum2 sum01 = 0.0f;

    auto [samples0, n0] = waveform0;
    auto [samples1, n1] = waveform1;

#ifdef MY_DEBUG
    if (n0 != n1) {
        printf("BUG 234f8273\n");
    }
#endif
    auto n = std::min(n0, n1);

    for (int64_t is = 0; is < n; ++is) {
        auto a0 = samples0[is];
        auto a1 = samples1[is];

        sum1 += a1;
        sum12 += a1*a1;
        sum01 += a0*a1;
    }

    {
        double nom = sum01*n - sum0*sum1;
        double den2a = sum02*n - sum0*sum0;
        double den2b = sum12*n - sum1*sum1;
        cc = (nom)/(sqrt(den2a*den2b));
    }

    return cc;
}

std::tuple<TCC, TOffset> findBestCC(
    const TWaveformView & waveform0,
    const TWaveformView & waveform1,
    int64_t alignWindow) {
    TCC bestcc = -1.0;
    TOffset besto = -1;

    auto [samples0, n0] = waveform0;
    auto [samples1, n1] = waveform1;

#ifdef MY_DEBUG
    if (n0 + 2*alignWindow != n1) {
        printf("BUG 924830jm92, n0 = %d, n1 = %d, a = %d\n", (int) n0, (int) n1, (int) alignWindow);
    }
#endif

    auto [sum0, sum02] = calcSum(waveform0);

    for (int o = 0; o < 2*alignWindow; ++o) {
        auto cc = calcCC(waveform0, { samples1 + o, n0 }, sum0, sum02);
        if (cc > bestcc) {
            besto = o - alignWindow;
            bestcc = cc;
        }
    }

    return { bestcc, besto };
}

bool calculateSimilartyMap(TKeyPressCollection & keyPresses, TSimilarityMap & res) {
    res.clear();
    int nPresses = keyPresses.size();

    int w = kKeyPressWidth_samples;
    int alignWindow = 256;

    res.resize(nPresses);
    for (auto & x : res) x.resize(nPresses);

    for (int i = 0; i < nPresses; ++i) {
        res[i][i].cc = 1.0f;
        res[i][i].offset = 0;

        auto & [waveform0, pos0, avgcc] = keyPresses[i];
        auto [samples0, n0] = waveform0;

        for (int j = 0; j < nPresses; ++j) {
            if (i == j) continue;

            auto waveform1 = keyPresses[j].waveform;
            auto pos1 = keyPresses[j].pos;

            auto samples1 = waveform1.samples;
            auto [bestcc, bestoffset] = findBestCC({ samples0 + pos0 + (int)(0.5f*w),               2*w },
                                                   { samples1 + pos1 + (int)(0.5f*w) - alignWindow, 2*w + 2*alignWindow }, alignWindow);

            res[i][i].cc = bestcc;
            res[i][i].offset = bestoffset;

            avgcc += bestcc;
        }
        avgcc /= (nPresses - 1);
    }

    return true;
}

void markKeyPresses(const TWaveform & waveform, const TKeyPressCollection & keyPresses, TWaveform & res) {
}

float plotWaveform(void * data, int i) {
    TWaveformView * waveform = (TWaveformView *)data;
    return waveform->samples[i];
};

int main(int argc, char ** argv) {
    srand(time(0));

    printf("Usage: %s record.kbd\n", argv[0]);
    if (argc < 2) {
        return -1;
    }

    if (SDL_Init(SDL_INIT_VIDEO|SDL_INIT_TIMER) != 0) {
        printf("Error: %s\n", SDL_GetError());
        return -1;
    }

    int windowSizeX = 600;
    int windowSizeY = 600;

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

    int64_t sampleRate = 24000;

    TWaveform waveformInput;
    printf("[+] Loading recording from '%s'\n", argv[1]);
    if (readFromFile(argv[1], waveformInput) == false) {
        printf("Specified file '%s' does not exist\n", argv[1]);
        return -1;
    }

    printf("[+] Loaded recording: of %d samples (sample size = %d bytes)\n", (int) waveformInput.size(), (int) sizeof(TSample));
    printf("    Size in memory:          %g MB\n", (float)(sizeof(TSample)*waveformInput.size())/1024/1024);
    printf("    Sample size:             %d\n", (int) sizeof(TSample));
    printf("    Total number of samples: %d\n", (int) waveformInput.size());
    printf("    Recording length:        %g seconds\n", (float)(waveformInput.size())/sampleRate);

    TKeyPressCollection keyPresses;
    {
        auto tStart = std::chrono::high_resolution_clock::now();
        printf("[+] Searching for key presses\n");
        if (findKeyPresses(waveformInput, keyPresses) == false) {
            printf("Failed to detect keypresses\n");
            return -2;
        }
        auto tEnd = std::chrono::high_resolution_clock::now();
        printf("[+] Detected a total of %d potential key presses\n", (int) keyPresses.size());
        for (auto & k : keyPresses) {
            printf("    position - %d\n", (int) k.pos);
        }
        printf("[+] Search took %4.3f seconds\n", toSeconds(tStart, tEnd));

        dumpKeyPresses("key_presses.plot", keyPresses);
    }

    TSimilarityMap similarityMap;
    {
        auto tStart = std::chrono::high_resolution_clock::now();
        printf("[+] Calculating CC similarity map\n");
        if (calculateSimilartyMap(keyPresses, similarityMap) == false) {
            printf("Failed to calculate similariy map\n");
            return -3;
        }
        auto tEnd = std::chrono::high_resolution_clock::now();
        printf("[+] Calculation took %4.3f seconds\n", toSeconds(tStart, tEnd));
    }

    int n = keyPresses.size();
    for (int i = 0; i < n; ++i) {
        auto avgcc = keyPresses[i].ccAvg;
        printf("    Average CC for keypress %4d - %6.3f\n", i, avgcc);
    }

    printf("%5d ", -1);
    for (int j = 0; j < n; ++j) {
        printf("%4d ", j);
    }
    printf("\n");
    printf("--------------------------------------------------------------------------------------------------------------------------------------\n");

    for (int i = 0; i < n; ++i) {
        printf("%2d  | ", i);
        for (int j = 0; j < n; ++j) {
            auto cc = similarityMap[i][j].cc;
            if (cc > -0.45) {
                printf("%4.0f ", cc*100);
            } else {
                printf("     ");
            }
            //printf("%3d - %3d -> %4.3f, %3d\n", i, j, cc, (int) offset);
        }
        printf("\n");
    }

    printf("\n");

    bool finishApp = false;
    while (finishApp == false) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            ImGui_ImplSDL2_ProcessEvent(&event);
            switch (event.type) {
                case SDL_QUIT:
                    finishApp = true;
                    break;
                case SDL_KEYDOWN:
                    if (event.key.keysym.sym == SDLK_ESCAPE) {
                        finishApp = true;
                    }
                    break;
                case SDL_WINDOWEVENT:
                    if (event.window.event == SDL_WINDOWEVENT_CLOSE && event.window.windowID == SDL_GetWindowID(window)) finishApp = true;
                    break;
            };
        }

        SDL_GetWindowSize(window, &windowSizeX, &windowSizeY);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplSDL2_NewFrame(window);
        ImGui::NewFrame();

        ImGui::Begin("Key Presses");
        static int nview = 1024;
        static int offset = 0;
        auto wview = getView(waveformInput, offset, nview);
        auto wsize = ImVec2(ImGui::GetContentRegionAvailWidth(), 250.0);
        auto savePos = ImGui::GetCursorScreenPos();
        auto drawList = ImGui::GetWindowDrawList();
        ImGui::PlotLines("##Waveform", plotWaveform, &wview, nview, 0, "Waveform", 16000, -16000, wsize);
        ImGui::SliderInt("View  ", &nview, 512, 48000);
        ImGui::DragInt("Offset", &offset, 0.01*nview, 0, waveformInput.size() - nview);
        offset = std::min((int) offset, (int) waveformInput.size() - nview);
        for (int i = 0; i < (int) keyPresses.size(); ++i) {
            if (keyPresses[i].pos + kKeyPressWidth_samples < offset) continue;
            if (keyPresses[i].pos - kKeyPressWidth_samples >= offset + nview) break;
            float x0 = ((float)(keyPresses[i].pos - kKeyPressWidth_samples - offset))/nview;
            float x1 = ((float)(keyPresses[i].pos + kKeyPressWidth_samples - offset))/nview;
            drawList->AddRectFilled({savePos.x + x0*wsize.x, savePos.y}, {savePos.x + x1*wsize.x, savePos.y + wsize.y}, ImGui::ColorConvertFloat4ToU32({1.0f, 0.0f, 0.0f, 0.3f}));
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

    return 0;
}
