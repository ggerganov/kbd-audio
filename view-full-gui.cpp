/*! \file view-full-gui.cpp
 *  \brief Visualize data recorded with record-full
 *  \author Georgi Gerganov
 */

#ifdef __EMSCRIPTEN__
#include "emscripten.h"
#define IMGUI_IMPL_OPENGL_LOADER_GLEW
#endif

#include "common.h"
#include "common-gui.h"
#include "constants.h"

#include "imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl3.h"

#include <SDL.h>

#if defined(IMGUI_IMPL_OPENGL_LOADER_GL3W)
#include <GL/gl3w.h>    // Initialize with gl3wInit()
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLEW)
#include <GL/glew.h>    // Initialize with glewInit()
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLAD)
#include <glad/glad.h>  // Initialize with gladLoadGL()
#else
#include IMGUI_IMPL_OPENGL_LOADER_CUSTOM
#endif

#include <chrono>
#include <cstdio>
#include <string>
#include <vector>
#include <thread>
#include <algorithm>
#include <functional>

static std::function<bool()> g_doInit;
static std::function<bool()> g_mainUpdate;

void mainUpdate() {
    g_mainUpdate();
}

// JS interface

extern "C" {
    int doInit() { return g_doInit(); }
}

#define MY_DEBUG

// globals
#ifdef __EMSCRIPTEN__
int g_windowSizeX = 1200;
int g_windowSizeY = 400;
#else
int g_windowSizeX = 1600;
int g_windowSizeY = 400;
#endif

struct stParameters;

using TParameters           = stParameters;

using TSampleInput          = TSampleF;
using TSample               = TSampleI16;
using TWaveform             = TWaveformI16;
using TWaveformView         = TWaveformViewI16;
using TPlaybackData         = TPlaybackDataI16;

struct stParameters {
    int playbackId              = 0;
    int keyPressWidth_samples   = 256;
    int sampleRate              = kSampleRate;
    int offsetFromPeak          = keyPressWidth_samples/2;
    int alignWindow             = 256;
    float thresholdClustering   = 0.5f;
};

float plotWaveform(void * data, int i) {
    TWaveformView * waveform = (TWaveformView *)data;
    return waveform->samples[i];
}

float plotWaveformInverse(void * data, int i) {
    TWaveformView * waveform = (TWaveformView *)data;
    return -waveform->samples[i];
}

SDL_AudioDeviceID g_deviceIdOut = 0;
TPlaybackData g_playbackData;

bool renderWaveform(TParameters & , const TWaveform & waveform) {
    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::SetNextWindowSize(ImVec2(g_windowSizeX, g_windowSizeY));
    if (ImGui::Begin("Waveform", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove)) {
        int viewMin = 512;
        int viewMax = waveform.size();

        static int nview = waveform.size();
        static int offset = (waveform.size() - nview)/2;
        static float amin = std::numeric_limits<TSample>::min();
        static float amax = std::numeric_limits<TSample>::max();
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
            ImGui::SliderFloat("Amplitude Min", &amin, std::numeric_limits<TSample>::min(), 0);
            ImGui::SliderFloat("Amplitude Max", &amax, 0, std::numeric_limits<TSample>::max());
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
                g_playbackData.playing = false;
                g_playbackData.idx = g_playbackData.waveform.n - TPlaybackData::kSamples;
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

        if (g_playbackData.idx >= g_playbackData.waveform.n) {
            g_playbackData.playing = false;
            SDL_ClearQueuedAudio(g_deviceIdOut);
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

    if (params.playbackId < 0 || params.playbackId >= nDevices) {
        printf("Invalid playback device id selected - %d\n", params.playbackId);
        return false;
    }

    SDL_AudioSpec playbackSpec;
    SDL_zero(playbackSpec);

    playbackSpec.freq = params.sampleRate;
    playbackSpec.format = std::is_same<TSample, int16_t>::value ? AUDIO_S16 : AUDIO_S32;
    playbackSpec.channels = 1;
    playbackSpec.samples = TPlaybackData::kSamples;
    playbackSpec.callback = cbPlayback<TSample>;
    playbackSpec.userdata = &g_playbackData;

    SDL_AudioSpec obtainedSpec;
    SDL_zero(obtainedSpec);

    g_deviceIdOut = SDL_OpenAudioDevice(SDL_GetAudioDeviceName(params.playbackId, SDL_FALSE), SDL_FALSE, &playbackSpec, &obtainedSpec, 0);
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
    printf("Usage: %s input.kbd [-pN]\n", argv[0]);
    printf("    -pN - select playback device N\n");
    printf("\n");

    printf("Usage: %s record.kbd\n", argv[0]);
    if (argc < 2) {
        return -1;
    }

    auto argm = parseCmdArguments(argc, argv);
    int playbackId = argm["p"].empty() ? 0 : std::stoi(argm["p"]);

    TParameters params;
    TWaveform waveformInput;

    if (SDL_Init(SDL_INIT_VIDEO|SDL_INIT_TIMER) != 0) {
        printf("Error: %s\n", SDL_GetError());
        return -1;
    }

    params.playbackId = playbackId;

    g_doInit = [&]() {
        return prepareAudioOut(params);
    };

    printf("[+] Loading recording from '%s'\n", argv[1]);
    if (readFromFile<TSampleF>(argv[1], waveformInput) == false) {
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
#ifdef __EMSCRIPTEN__
	SDL_Window* window;
	SDL_Renderer *renderer;
    SDL_CreateWindowAndRenderer(g_windowSizeX, g_windowSizeY, SDL_WINDOW_OPENGL, &window, &renderer);
#else
    SDL_Window* window = SDL_CreateWindow("View-full", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, g_windowSizeX, g_windowSizeY, SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE|SDL_WINDOW_ALLOW_HIGHDPI);
#endif
    SDL_GLContext gl_context = SDL_GL_CreateContext(window);
    SDL_GL_SetSwapInterval(1); // Enable vsync

    // Initialize OpenGL loader
#if defined(IMGUI_IMPL_OPENGL_LOADER_GL3W)
    bool err = gl3wInit() != 0;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLEW)
    bool err = glewInit() != GLEW_OK;
#elif defined(IMGUI_IMPL_OPENGL_LOADER_GLAD)
    bool err = gladLoadGL() == 0;
#else
    bool err = false; // If you use IMGUI_IMPL_OPENGL_LOADER_CUSTOM, your loader is likely to requires some form of initialization.
#endif
    if (err)
    {
        fprintf(stderr, "Failed to initialize OpenGL loader!\n");
        return 1;
    }

    // Setup Dear ImGui binding
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard Controls

    ImGui_ImplSDL2_InitForOpenGL(window, gl_context);
#ifdef __EMSCRIPTEN__
    ImGui_ImplOpenGL3_Init("#version 300 es");
#else
    ImGui_ImplOpenGL3_Init(glsl_version);
#endif

    // Setup style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsClassic();

    ImFontConfig fontConfig;
    //fontConfig.SizePixels = 14.0f;
    ImGui::GetIO().Fonts->AddFontDefault(&fontConfig);

    ImGui::GetStyle().AntiAliasedFill = false;
    ImGui::GetStyle().AntiAliasedLines = false;

    printf("[+] Loaded recording: of %d samples (sample size = %d bytes)\n", (int) waveformInput.size(), (int) sizeof(TSample));
    printf("    Size in memory:          %g MB\n", (float)(sizeof(TSample)*waveformInput.size())/1024/1024);
    printf("    Sample size:             %d\n", (int) sizeof(TSample));
    printf("    Total number of samples: %d\n", (int) waveformInput.size());
    printf("    Recording length:        %g seconds\n", (float)(waveformInput.size())/params.sampleRate);

    bool finishApp = false;
    g_mainUpdate = [&]() {
        if (finishApp) return false;
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

        renderWaveform(params, waveformInput);

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

        return true;
    };

#ifdef __EMSCRIPTEN__
    emscripten_set_main_loop(mainUpdate, 60, 1);
#else
    if (g_doInit() == false) {
        printf("Error: failed to initialize audio playback\n");
        return -2;
    }

    while (true) {
        if (g_mainUpdate() == false) break;
    }
#endif

    printf("[+] Terminated");

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplSDL2_Shutdown();
    ImGui::DestroyContext();

    SDL_GL_DeleteContext(gl_context);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
