/*! \file key_average_gui.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "key_logger.h"
#include "audio_logger.h"

#include "fftw3.h"

#include "imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl3.h"

#include <SDL.h>
#include <GL/gl3w.h>

#include <map>
#include <mutex>
#include <string>

int main(int, char**) {
    if (SDL_Init(SDL_INIT_VIDEO|SDL_INIT_AUDIO|SDL_INIT_TIMER) != 0) {
        printf("Error: %s\n", SDL_GetError());
        return -1;
    }

    int windowSizeX = 1920;
    int windowSizeY = 1080;

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
    SDL_Window* window = SDL_CreateWindow("Average key stroke audio information", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, windowSizeX, windowSizeY, SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE);
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

    // Project specific
    std::mutex mutex;

    auto fftIn = (float*) fftwf_malloc(sizeof(float)*AudioLogger::kSamplesPerFrame);
    auto fftOut = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*AudioLogger::kSamplesPerFrame);
    auto fftPlan = fftwf_plan_dft_r2c_1d(1*AudioLogger::kSamplesPerFrame, fftIn, fftOut, FFTW_ESTIMATE);

    using TKey = char;
    using TKeyAverage = std::array<AudioLogger::Frame, 2*AudioLogger::kBufferSize_frames - 1>;

    TKey keyPressed = -1;
    std::map<TKey, TKeyAverage> keySoundAverageAmpl;
    std::map<TKey, TKeyAverage> keySoundAverageFreq;

    AudioLogger audioLogger;
    AudioLogger::Callback cbAudio = [&](const auto & frames) {
        std::lock_guard<std::mutex> lock(mutex);

        int fid = 0;
        auto & buffersAmpl = keySoundAverageAmpl[keyPressed];
        auto & buffersFreq = keySoundAverageFreq[keyPressed];
        for (const auto & frame : frames) {
            for (auto i = 0; i < AudioLogger::kSamplesPerFrame; ++i) {
                buffersAmpl[fid][i] += frame[i];
            }
            std::copy(buffersAmpl[fid].begin(), buffersAmpl[fid].end(), fftIn);
            fftwf_execute(fftPlan);
            for (auto i = 0; i < AudioLogger::kSamplesPerFrame; ++i) {
                buffersFreq[fid][i] = fftOut[i][0]*fftOut[i][0] + fftOut[i][1]*fftOut[i][1];
            }
            for (auto i = 0; i < AudioLogger::kSamplesPerFrame/2; ++i) {
                buffersFreq[fid][i] += buffersFreq[fid][AudioLogger::kSamplesPerFrame - i];
            }

            if (++fid >= buffersAmpl.size()) break;
        }

        keyPressed = -1;
    };

    if (audioLogger.install(cbAudio) == false) {
        fprintf(stderr, "Failed to install audio logger\n");
        return -1;
    }

    KeyLogger keyLogger;
    KeyLogger::Callback cbKey = [&audioLogger, &keyPressed](int key) -> void {
        if (keyPressed == -1) {
            keyPressed = key;
            audioLogger.record();
        }
    };

    if (keyLogger.install(cbKey) == false) {
        fprintf(stderr, "Failed to install key logger\n");
        return -2;
    }

    // Main loop
    bool done = false;
    while (!done) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            ImGui_ImplSDL2_ProcessEvent(&event);
            if (event.type == SDL_QUIT)
                done = true;
            if (event.type == SDL_WINDOWEVENT && event.window.event == SDL_WINDOWEVENT_CLOSE && event.window.windowID == SDL_GetWindowID(window))
                done = true;
        }

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplSDL2_NewFrame(window);
        ImGui::NewFrame();

        ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Once);
        ImGui::SetNextWindowSize(ImVec2(windowSizeX, windowSizeY), ImGuiCond_Once);
        ImGui::Begin("Average Key Waveform");
        ImGui::Text("Frames in buffer: %d\n", (int) AudioLogger::kBufferSize_frames);
        {
            std::lock_guard<std::mutex> lock(mutex);

            for (const auto & ampl : keySoundAverageAmpl) {
                struct SampleGetter {
                    static float f(void * data, int i) {
                        int fid = i/AudioLogger::kSamplesPerFrame;
                        int sid = i - fid*AudioLogger::kSamplesPerFrame;
                        const auto & frames = *(TKeyAverage *)(data);
                        const auto & frame = frames[fid];
                        return frame[sid];
                    }
                };

                float (*getter)(void *, int) = SampleGetter::f;

                std::string skey(KeyLogger::codeToText(ampl.first));
                ImGui::PlotLines(
                    skey.c_str(),
                    getter, (void *)(intptr_t)(&ampl.second), ampl.second.size()*AudioLogger::kSamplesPerFrame,
                    0, skey.c_str(), FLT_MAX, FLT_MAX, ImVec2(windowSizeX, 0.1f*windowSizeY));

                int nFrames = keySoundAverageFreq[ampl.first].size();
                for (int fid = 0; fid < (int) keySoundAverageFreq[ampl.first].size(); ++fid) {
                    const auto & freq = keySoundAverageFreq[ampl.first][fid];
                    ImGui::PlotHistogram(
                        ("##" + skey + " frame " + std::to_string(fid)).c_str(),
                        freq.data(), freq.size()/2,
                        0, ("##" + skey + " frame " + std::to_string(fid)).c_str(), 0.0f, 50.0f, ImVec2(windowSizeX/nFrames, 0.1f*windowSizeY));
                    ImGui::SameLine();
                }
                ImGui::Text("%s", "");
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

    // Cleanup
    if (fftPlan) fftwf_destroy_plan(fftPlan);
    if (fftIn) fftwf_free(fftIn);
    if (fftOut) fftwf_free(fftOut);

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplSDL2_Shutdown();
    ImGui::DestroyContext();

    SDL_GL_DeleteContext(gl_context);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
