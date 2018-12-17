/*! \file view-gui.cpp
 *  \brief Visualize training data recorded with 'record'
 *  \author Georgi Gerganov
 */

#include "constants.h"

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
#include <thread>
#include <algorithm>

#define MY_DEBUG

// globals
int g_windowSizeX = 1600;
int g_windowSizeY = 400;

struct stParameters;
struct stWaveformView;

using TKey                  = int32_t;
using TParameters           = stParameters;

using TSampleInput          = float;
using TSample               = int32_t;
using TTrainKeys            = std::vector<TKey>;
using TWaveform             = std::vector<TSample>;
using TWaveformView         = stWaveformView;

struct stParameters {
    int keyPressWidth_samples   = 256;
    int sampleRate              = kSampleRate;
    int offsetFromPeak          = keyPressWidth_samples/2;
    int alignWindow             = 256;
    float thresholdClustering   = 0.5f;
};

struct stWaveformView {
    const TSample * samples     = nullptr;
    int64_t                     n = 0;
};

template <typename T>
float toSeconds(T t0, T t1) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()/1024.0f;
}

TWaveformView getView(const TWaveform & waveform, int64_t idx) { return { waveform.data() + idx, (int64_t) waveform.size() - idx }; }

TWaveformView getView(const TWaveform & waveform, int64_t idx, int64_t len) { return { waveform.data() + idx, len }; }

bool readFromFile(const TParameters & params, const std::string & fname, TWaveform & res, TTrainKeys & trainKeys) {
    trainKeys.clear();

    std::ifstream fin(fname, std::ios::binary);
    if (fin.good() == false) {
        return false;
    }

    int32_t bufferSize_frames = 1;
    fin.read((char *)(&bufferSize_frames), sizeof(bufferSize_frames));
    if (bufferSize_frames != kTrainBufferSize_frames) {
        printf("Buffer size in file (%d) does not match the expected one (%d)\n", bufferSize_frames, (int) kTrainBufferSize_frames);
        return false;
    }

    {
        static_assert(std::is_same<TSampleInput, float>::value, "TSampleInput not recognised");
        static_assert(
            std::is_same<TSample, float>::value
            || std::is_same<TSample, int16_t>::value
            || std::is_same<TSample, int32_t>::value
                      , "TSampleInput not recognised");

        int32_t offset = 0;
        std::streamsize size = bufferSize_frames*kSamplesPerFrame*sizeof(TSampleInput);
        while (true) {
            TKey keyPressed = 0;
            fin.read((char *)(&keyPressed), sizeof(keyPressed));
            if (fin.eof()) break;
            trainKeys.push_back(keyPressed);

            if (std::is_same<TSample, int16_t>::value) {
                std::vector<TSampleInput> buf(size/sizeof(TSampleInput));
                res.resize(offset + size/sizeof(TSampleInput));
                fin.read((char *)(buf.data()), size);
                double amax = 0.0f;
                for (auto i = 0; i < buf.size(); ++i) if (std::abs(buf[i]) > amax) amax = std::abs(buf[i]);
                for (auto i = 0; i < buf.size(); ++i) res[offset + i] = std::round(32000.0*(buf[i]/amax));
            } else if (std::is_same<TSample, int32_t>::value) {
                std::vector<TSampleInput> buf(size/sizeof(TSampleInput));
                res.resize(offset + size/sizeof(TSampleInput));
                fin.read((char *)(buf.data()), size);
                double amax = 0.0f;
                for (auto i = 0; i < buf.size(); ++i) if (std::abs(buf[i]) > amax) amax = std::abs(buf[i]);
                for (auto i = 0; i < buf.size(); ++i) res[offset + i] = std::round(32000.0*(buf[i]/amax));
            } else if (std::is_same<TSample, float>::value) {
                res.resize(offset + size/sizeof(TSample));
                fin.read((char *)(res.data() + offset), size);
            } else {
            }

            offset += size/sizeof(TSampleInput);
            if (fin.eof()) break;
        }
    }

    fin.close();

    return true;
}

bool generateLowResWaveform(const TWaveformView & waveform, TWaveform & waveformLowRes, int nWindow) {
    waveformLowRes.resize(waveform.n);

    int k = nWindow;
    std::deque<int64_t> que(k);

    //auto [samples, n] = waveform;
    auto samples = waveform.samples;
    auto n       = waveform.n;

    TWaveform waveformAbs(n);
    for (int64_t i = 0; i < n; ++i) {
        waveformAbs[i] = std::abs(samples[i]);
    }

    for (int64_t i = 0; i < n; ++i) {
        if (i < k) {
            while((!que.empty()) && waveformAbs[i] >= waveformAbs[que.back()]) {
                que.pop_back();
            }
            que.push_back(i);
        } else {
            while((!que.empty()) && que.front() <= i - k) {
                que.pop_front();
            }

            while((!que.empty()) && waveformAbs[i] >= waveformAbs[que.back()]) {
                que.pop_back();
            }

            que.push_back(i);

            int64_t itest = i - k/2;
            waveformLowRes[itest] = waveformAbs[que.front()];
        }
    }

    return true;
}

bool generateLowResWaveform(const TWaveform & waveform, TWaveform & waveformLowRes, int nWindow) {
    return generateLowResWaveform(getView(waveform, 0), waveformLowRes, nWindow);
}

float plotWaveform(void * data, int i) {
    TWaveformView * waveform = (TWaveformView *)data;
    return waveform->samples[i];
};

float plotWaveformInverse(void * data, int i) {
    TWaveformView * waveform = (TWaveformView *)data;
    return -waveform->samples[i];
};

struct PlaybackData {
    static const int kSamples = 1024;
    bool playing = false;
    int slowDown = 1;
    int64_t idx = 0;
    int64_t offset = 0;
    TWaveformView waveform;
};

SDL_AudioDeviceID g_deviceIdOut = 0;
PlaybackData g_playbackData;

bool renderWaveform(TParameters & params, const TWaveform & waveform, const TTrainKeys & trainKeys) {
    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::SetNextWindowSize(ImVec2(g_windowSizeX, g_windowSizeY));
    if (ImGui::Begin("Waveform", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove)) {
        int viewMin = 512;
        int viewMax = waveform.size();

        static int nview = waveform.size();
        static int offset = (waveform.size() - nview)/2;
        static float amin = -16000;
        static float amax = 16000;
        static float dragOffset = 0.0f;
        static float scrollSize = 18.0f;

        static auto nviewPrev = nview + 1;

        static TWaveform waveformLowRes = waveform;
        static TWaveform waveformThreshold = waveform;

        auto wview = getView(waveformLowRes, offset, nview);
        auto wsize = ImGui::GetContentRegionAvail();
        wsize.y -= 50.0f;

        auto mpos = ImGui::GetIO().MousePos;
        auto savePos = ImGui::GetCursorScreenPos();
        auto drawList = ImGui::GetWindowDrawList();
        ImGui::PushStyleColor(ImGuiCol_FrameBg, { 0.3f, 0.3f, 0.3f, 0.3f });
        ImGui::PushStyleColor(ImGuiCol_PlotHistogram, { 1.0f, 1.0f, 1.0f, 1.0f });
        ImGui::PlotHistogram("##Waveform", plotWaveformInverse, &wview, nview, 0, "Waveform", amin, amax, wsize);
        ImGui::PopStyleColor(2);
        ImGui::SetCursorScreenPos(savePos);
        ImGui::PushStyleColor(ImGuiCol_FrameBg, { 0.1f, 0.1f, 0.1f, 0.0f });
        ImGui::PushStyleColor(ImGuiCol_PlotHistogram, { 1.0f, 1.0f, 1.0f, 1.0f });
        ImGui::PlotHistogram("##Waveform", plotWaveform, &wview, nview, 0, "Waveform", amin, amax, wsize);
        ImGui::PopStyleColor(2);
        ImGui::SetCursorScreenPos(savePos);
        ImGui::InvisibleButton("##WaveformIB",wsize);
        if (ImGui::IsItemHovered()) {
            auto w = ImGui::GetIO().MouseWheel;
            auto nview_old = nview;
            nview *= (10.0 + w)/10.0;
            nview = std::min(std::max(viewMin, nview), viewMax);
            if (w != 0.0) {
                offset = std::max(0.0f, offset + ((mpos.x - savePos.x)/wsize.x)*(nview_old - nview));
            }

            if (ImGui::IsMouseDown(0) && ImGui::IsMouseDragging(0) == false) {
                dragOffset = offset;
            }

            if (ImGui::IsMouseDragging(0)) {
                offset = dragOffset - ImGui::GetMouseDragDelta(0).x*nview/wsize.x;
            }
        }
        if (ImGui::BeginPopupContextWindow()) {
            ImGui::SliderInt("View  ", &nview, viewMin, viewMax);
            ImGui::DragInt  ("Offset", &offset, 0.01*nview, 0, waveform.size() - nview);
            ImGui::SliderFloat("Amplitude Min", &amin, -32000, 0);
            ImGui::SliderFloat("Amplitude Max", &amax, 0, 32000);
            ImGui::EndPopup();
        }

        ImGui::InvisibleButton("##WaveformScrollIB", {wsize.x, scrollSize});
        drawList->AddRect({savePos.x, savePos.y + wsize.y}, {savePos.x + wsize.x, savePos.y + wsize.y + scrollSize}, ImGui::ColorConvertFloat4ToU32({1.0f, 1.0f, 1.0f, 1.0f}));
        drawList->AddRectFilled({savePos.x + wsize.x*(1.f*offset)/viewMax, savePos.y + wsize.y}, {savePos.x + wsize.x*(1.f*offset + nview)/viewMax, savePos.y + wsize.y + scrollSize}, ImGui::ColorConvertFloat4ToU32({1.0f, 1.0f, 1.0f, 1.0f}));

        auto savePos2 = ImGui::GetCursorScreenPos();

        static bool scrolling = false;
        if (ImGui::IsItemHovered()) {
            if (ImGui::IsMouseDown(0)) {
                scrolling = true;
            }
        }

        if (scrolling) {
            offset = ((mpos.x - savePos.x)/wsize.x)*viewMax - nview/2;
        }

        if (ImGui::IsMouseDown(0) == false) {
            scrolling = false;
        }

        offset = std::max(0, std::min((int) offset, (int) waveform.size() - nview));

        for (int i = 0; i < (int) trainKeys.size(); ++i) {
            int pos = i*kSamplesPerFrame*5 + 2.5*kSamplesPerFrame;
            if (pos + params.offsetFromPeak + params.keyPressWidth_samples < offset) continue;
            if (pos + params.offsetFromPeak - params.keyPressWidth_samples >= offset + nview) break;

            {
                float x0 = ((float)(pos - offset))/nview;

                ImVec2 p0 = { savePos.x + x0*wsize.x, savePos.y };
                ImVec2 p1 = { savePos.x + x0*wsize.x, savePos.y + wsize.y };

                drawList->AddLine(p0, p1, ImGui::ColorConvertFloat4ToU32({ 1.0f, 0.0f, 0.0f, 0.8f }), 1.0f);
            }

            {
                float x0 = ((float)(pos - offset))/nview;
                float x1 = ((float)(pos + params.offsetFromPeak - params.keyPressWidth_samples - offset))/nview;
                float x2 = ((float)(pos + params.offsetFromPeak + params.keyPressWidth_samples - offset))/nview;

                ImVec2 p0 = { savePos.x + x0*wsize.x, savePos.y };
                ImVec2 p1 = { savePos.x + x1*wsize.x, savePos.y };
                ImVec2 p2 = { savePos.x + x2*wsize.x, savePos.y + wsize.y };

                if (nview < 64.0*wsize.x) {
                    ImGui::SetCursorScreenPos({ savePos.x + 0.5f*((x1 + x2)*wsize.x - ImGui::CalcTextSize(std::to_string(i).c_str()).x), savePos.y + wsize.y - ImGui::GetTextLineHeightWithSpacing() });
                    ImGui::Text("%d - %s", i, kKeyText.at(trainKeys[i]));
                }

            }
        }

        {
            float x0 = ((float)(g_playbackData.offset + g_playbackData.idx - offset))/nview;

            ImVec2 p0 = {savePos.x + x0*wsize.x, savePos.y};
            ImVec2 p1 = {savePos.x + x0*wsize.x, savePos.y + wsize.y};
            drawList->AddLine(p0, p1, ImGui::ColorConvertFloat4ToU32({ 1.0f, 1.0f, 0.0f, 0.3f }));
        }

        ImGui::SetCursorScreenPos(savePos2);

        //auto io = ImGui::GetIO();
        //ImGui::Text("Keys pressed:");   for (int i = 0; i < IM_ARRAYSIZE(io.KeysDown); i++) if (ImGui::IsKeyPressed(i))             { ImGui::SameLine(); ImGui::Text("%d", i); }

        static bool playHalfSpeed = false;
        ImGui::PushItemWidth(100.0);

        ImGui::Checkbox("x0.5", &playHalfSpeed);
        ImGui::SameLine();
        if (g_playbackData.playing) {
            if (ImGui::Button("Stop") || ImGui::IsKeyPressed(44)) { // space
                g_playbackData.idx = g_playbackData.waveform.n - 1;
            }
        } else {
            if (ImGui::Button("Play") || ImGui::IsKeyPressed(44)) { // space
                g_playbackData.playing = true;
                g_playbackData.slowDown = playHalfSpeed ? 2 : 1;
                g_playbackData.idx = 0;
                g_playbackData.offset = offset;
                g_playbackData.waveform = getView(waveform, offset, nview);
                SDL_PauseAudioDevice(g_deviceIdOut, 0);
            }
        }

        if (g_playbackData.idx > g_playbackData.waveform.n) {
            g_playbackData.playing = false;
            SDL_PauseAudioDevice(g_deviceIdOut, 1);
        }

        ImGui::PopItemWidth();

        if (nview != nviewPrev) {
            generateLowResWaveform(waveform, waveformLowRes, std::max(1.0f, nview/wsize.x));
            nviewPrev = nview;
        }
    }
    ImGui::End();

    return false;
}

void cbPlayback(void * userData, uint8_t * stream, int len) {
    PlaybackData * data = (PlaybackData *)(userData);
    auto end = std::min(data->idx + PlaybackData::kSamples/data->slowDown, data->waveform.n);
    auto idx = data->idx;
    auto sidx = 0;
    for (; idx < end; ++idx) {
        int16_t a = data->waveform.samples[idx];
        memcpy(stream + (sidx)*sizeof(a), &a, sizeof(a));
        len -= sizeof(a);
        ++sidx;

        if (data->slowDown == 2) {
            int16_t a2 = data->waveform.samples[idx + 1];
            a = 0.5*(a + a2);
            memcpy(stream + (sidx)*sizeof(a), &a, sizeof(a));
            len -= sizeof(a);
            ++sidx;
        }
    }
    while (len > 0) {
        int16_t a = 0;
        memcpy(stream + (idx - data->idx)*sizeof(a), &a, sizeof(a));
        len -= sizeof(a);
        ++idx;
    }
    data->idx = idx;
}

bool prepareAudioOut(const TParameters & params) {
    if (SDL_Init(SDL_INIT_AUDIO) < 0) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't initialize SDL: %s\n", SDL_GetError());
        return false;
    }

    int nDevices = SDL_GetNumAudioDevices(SDL_FALSE);
    printf("Found %d playback devices:\n", nDevices);
    for (int i = 0; i < nDevices; i++) {
        printf("    - Playback device #%d: '%s'\n", i, SDL_GetAudioDeviceName(i, SDL_FALSE));
    }

    SDL_AudioSpec playbackSpec;
    SDL_zero(playbackSpec);

    playbackSpec.freq = params.sampleRate;
    playbackSpec.format = AUDIO_S16;
    playbackSpec.channels = 1;
    playbackSpec.samples = PlaybackData::kSamples;
    playbackSpec.callback = cbPlayback;
    playbackSpec.userdata = &g_playbackData;

    SDL_AudioSpec obtainedSpec;
    SDL_zero(obtainedSpec);

    g_deviceIdOut = SDL_OpenAudioDevice(NULL, SDL_FALSE, &playbackSpec, &obtainedSpec, 0);
    if (!g_deviceIdOut) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't open an audio device for playback: %s!\n", SDL_GetError());
        SDL_Quit();
        return false;
    }

    printf("Opened playback device succesfully!\n");
    printf("    Frequency:  %d\n", obtainedSpec.freq);
    printf("    Format:     %d\n", obtainedSpec.format);
    printf("    Channels:   %d\n", obtainedSpec.channels);
    printf("    Samples:    %d\n", obtainedSpec.samples);

    SDL_PauseAudioDevice(g_deviceIdOut, 1);

    return true;
}

int main(int argc, char ** argv) {
    srand(time(0));

    printf("Usage: %s record.kbd\n", argv[0]);
    if (argc < 2) {
        return -1;
    }

    TParameters params;
    TWaveform waveformInput;
    TTrainKeys trainKeys;

    if (SDL_Init(SDL_INIT_VIDEO|SDL_INIT_TIMER) != 0) {
        printf("Error: %s\n", SDL_GetError());
        return -1;
    }

    if (prepareAudioOut(params) == false) {
        printf("Error: failed to initialize audio playback\n");
        return -2;
    }

    printf("[+] Loading recording from '%s'\n", argv[1]);
    if (readFromFile(params, argv[1], waveformInput, trainKeys) == false) {
        printf("Specified file '%s' does not exist\n", argv[1]);
        return -1;
    }

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
    SDL_Window* window = SDL_CreateWindow("View", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, g_windowSizeX, g_windowSizeY, SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE|SDL_WINDOW_ALLOW_HIGHDPI);
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

    ImGui::GetStyle().AntiAliasedFill = false;
    ImGui::GetStyle().AntiAliasedLines = false;

    ImGui_ImplSDL2_InitForOpenGL(window, gl_context);
    ImGui_ImplOpenGL3_Init(glsl_version);

    printf("[+] Loaded recording: of %d samples (sample size = %d bytes)\n", (int) waveformInput.size(), (int) sizeof(TSample));
    printf("    Size in memory:          %g MB\n", (float)(sizeof(TSample)*waveformInput.size())/1024/1024);
    printf("    Sample size:             %d\n", (int) sizeof(TSample));
    printf("    Total number of samples: %d\n", (int) waveformInput.size());
    printf("    Recording length:        %g seconds\n", (float)(waveformInput.size())/params.sampleRate);

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

        SDL_GetWindowSize(window, &g_windowSizeX, &g_windowSizeY);

        auto tStart = std::chrono::high_resolution_clock::now();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplSDL2_NewFrame(window);
        ImGui::NewFrame();

        renderWaveform(params, waveformInput, trainKeys);

        ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);

        ImGui::Render();
        SDL_GL_MakeCurrent(window, gl_context);
        glViewport(0, 0, (int) io.DisplaySize.x, (int) io.DisplaySize.y);
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        SDL_GL_SwapWindow(window);

        // stupid hack to limit frame-rate to ~60 fps on Mojave
        auto tEnd = std::chrono::high_resolution_clock::now();
        auto tus = std::chrono::duration_cast<std::chrono::microseconds>(tEnd - tStart).count();
        while (tus < 1e6/60.0) {
            std::this_thread::sleep_for(std::chrono::microseconds(std::max(100, (int) (0.5*(1e6/60.0 - tus)))));
            tEnd = std::chrono::high_resolution_clock::now();
            tus = std::chrono::duration_cast<std::chrono::microseconds>(tEnd - tStart).count();
        }
    }

    return 0;
}
