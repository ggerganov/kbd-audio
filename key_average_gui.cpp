/*! \file key_average_gui.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "constants.h"
#include "common.h"
#include "audio_logger.h"

#include "fftw3.h"

#include "imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl3.h"

#include <SDL.h>
#include <GL/gl3w.h>

#include <map>
#include <mutex>
#include <cmath>
#include <string>
#include <cstdlib>
#include <chrono>

int main(int argc, char ** argv) {
    if (SDL_Init(SDL_INIT_VIDEO|SDL_INIT_AUDIO|SDL_INIT_TIMER) != 0) {
        printf("Error: %s\n", SDL_GetError());
        return -1;
    }

    auto argm = parseCmdArguments(argc, argv);
    int captureId = argm["c"].empty() ? 0 : std::stoi(argm["c"]);
    int nChannels = argm["C"].empty() ? 0 : std::stoi(argm["C"]);

    printf("Usage: %s output.kbd [-cN]\n", argv[0]);
    printf("    -cN - select capture device N\n");
    printf("    -CN - number N of capture channels N\n");
    printf("\n");

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

    printf("sizeof fftwf_complex = %d\n", (int) sizeof(fftwf_complex));
    fftwf_complex * fftIn = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*kSamplesPerFrame);
    fftwf_complex * fftOut = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*kSamplesPerFrame);
    fftwf_plan fftPlan = fftwf_plan_dft_1d(kSamplesPerFrame, fftIn, fftOut, FFTW_FORWARD, FFTW_ESTIMATE);

    using TKey = int;
    using TKeyWaveform = std::array<AudioLogger::Frame, kBufferSizeTrain_frames>;

    TKey keyPressed = -1;
    std::map<TKey, TKeyWaveform> keySoundAverageAmpl;
    std::map<TKey, TKeyWaveform> keySoundAverageFreq;

    std::vector<float> similarityForOffset(4*kSamplesPerFrame);

    AudioLogger audioLogger;
    AudioLogger::Callback cbAudio = [&](const AudioLogger::Record & framesOriginal) {
        //auto t1 = std::chrono::high_resolution_clock::now();
        //printf("Received: %d\n", (int) std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());

        auto frames = framesOriginal;

        int nFrames = frames.size();
        bool exist = keySoundAverageAmpl.find(keyPressed) != keySoundAverageAmpl.end();
        exist = false;
        if (exist) {
            int kmax = 0;
            int imax = 0;
            float amax = 0.0f;

            int n = kSamplesPerFrame;
            for (int k = 0; k < nFrames; ++k) {
                auto & curFrame = frames[k];
                for (int i = 0; i < n; ++i) {
                    float a = std::abs(curFrame[i]);
                    if (a > amax) {
                        kmax = k;
                        imax = i;
                        amax = a;
                    }
                }
            }

            int besto = 0;
            float bestcc = 0.0f;
            float bestccq = 0.0f;
            float bestccp = 0.0f;

            std::fill(similarityForOffset.begin(), similarityForOffset.end(), 0);

            for (int o = -n; o < n; ++o) {
                int w = n;

                int s00 = kmax*n + imax - w;
                int s01 = kmax*n + imax + w;
                if (s00 < 0 || s01 >= nFrames*n) continue;

                int s10 = kmax*n + imax + o - w;
                int s11 = kmax*n + imax + o + w;
                if (s10 < 0 || s11 >= nFrames*n) continue;

                float cc = 0.0f;
                float ccq = 0.0f;
                float ccp = 0.0f;

                float suma = 0.0f, suma2 = 0.0f, sumb = 0.0f, sumb2 = 0.0f, sumab = 0.0f;
                float sumaq = 0.0f, suma2q = 0.0f, sumabq = 0.0f;
                float sumap = 0.0f, suma2p = 0.0f, sumabp = 0.0f;
                for (int s0 = s00; s0 < s01; ++s0) {
                    int k0 = s0/n;
                    int i0 = s0%n;

                    int s1 = s0 + o;
                    int k1 = s1/n;
                    int i1 = s1%n;

                    sumb += frames[k1][i1];
                    sumb2 += frames[k1][i1]*frames[k1][i1];

                    suma += keySoundAverageAmpl[keyPressed][k0][i0];
                    suma2 += keySoundAverageAmpl[keyPressed][k0][i0]*keySoundAverageAmpl[keyPressed][k0][i0];
                    sumab += frames[k1][i1]*keySoundAverageAmpl[keyPressed][k0][i0];
                    sumaq += keySoundAverageAmpl['q'][k0][i0];
                    suma2q += keySoundAverageAmpl['q'][k0][i0]*keySoundAverageAmpl['q'][k0][i0];
                    sumabq += frames[k1][i1]*keySoundAverageAmpl['q'][k0][i0];
                    sumap += keySoundAverageAmpl['p'][k0][i0];
                    suma2p += keySoundAverageAmpl['p'][k0][i0]*keySoundAverageAmpl['p'][k0][i0];
                    sumabp += frames[k1][i1]*keySoundAverageAmpl['p'][k0][i0];

                    //cc += frames[k1][i1]*keySoundAverageAmpl[keyPressed][k0][i0];
                    //ccq += frames[k1][i1]*keySoundAverageAmpl['q'][k0][i0];
                    //ccp += frames[k1][i1]*keySoundAverageAmpl['p'][k0][i0];
                }

                int ncc = s01 - s00;
                {
                    float nom = sumab*ncc - suma*sumb;
                    float den2a = suma2*ncc - suma*suma;
                    float den2b = sumb2*ncc - sumb*sumb;
                    cc = (nom)/(sqrt(den2a*den2b));
                }

                {
                    float nom = sumabq*ncc - sumaq*sumb;
                    float den2a = suma2q*ncc - sumaq*sumaq;
                    float den2b = sumb2*ncc - sumb*sumb;
                    ccq = (nom)/(sqrt(den2a*den2b));
                }

                {
                    float nom = sumab*ncc - sumap*sumb;
                    float den2a = suma2p*ncc - sumap*sumap;
                    float den2b = sumb2*ncc - sumb*sumb;
                    ccp = (nom)/(sqrt(den2a*den2b));
                }

                similarityForOffset[o + 2*n] = cc;

                if (cc > bestcc) {
                    besto = o;
                    bestcc = cc;
                }

                if (ccq > bestccq) {
                    bestccq = ccq;
                }

                if (ccp > bestccp) {
                    bestccp = ccp;
                }
            }

            printf("besto = %d\n", besto);
            printf("more similar to '%c'\n", (bestccq > bestccp) ? 'q' : 'p');

            if (besto > 0) {
                for (int s0 = 0; s0 < nFrames*n; ++s0) {
                    int k0 = s0/n;
                    int i0 = s0%n;

                    int s1 = s0 + besto;
                    int k1 = s1/n;
                    int i1 = s1%n;
                    frames[k0][i0] = s1 < nFrames*n ? frames[k1][i1] : 0.0f;
                }
            } else {
                for (int s0 = nFrames*n - 1; s0 >= 0; --s0) {
                    int k0 = s0/n;
                    int i0 = s0%n;

                    int s1 = s0 + besto;
                    int k1 = s1/n;
                    int i1 = s1%n;
                    frames[k0][i0] = s1 >= 0 ? frames[k1][i1] : 0.0f;
                }
            }
        }

        std::lock_guard<std::mutex> lock(mutex);

        int fid = 0;
        auto & buffersAmpl = keySoundAverageAmpl[keyPressed];
        auto & buffersFreq = keySoundAverageFreq[keyPressed];
        for (int k = 0; k < nFrames; ++k) {
            for (auto i = 0; i < kSamplesPerFrame; ++i) {
                buffersAmpl[fid][i] = frames[k][i];
                fftIn[i][0] = buffersAmpl[fid][i];
                fftIn[i][1] = 0;
            }
            fftwf_execute(fftPlan);
            for (auto i = 0; i < kSamplesPerFrame; ++i) {
                buffersFreq[fid][i] = fftOut[i][0]*fftOut[i][0] + fftOut[i][1]*fftOut[i][1];
            }
            for (auto i = 1; i < kSamplesPerFrame/2; ++i) {
                buffersFreq[fid][i] += buffersFreq[fid][kSamplesPerFrame - i];
            }

            if (++fid >= (int) buffersAmpl.size()) break;
        }

        keyPressed = -1;
    };

    AudioLogger::Parameters parameters;
    parameters.callback = std::move(cbAudio);
    parameters.captureId = captureId;
    parameters.nChannels = nChannels;
    parameters.sampleRate = kSampleRate;
    parameters.freqCutoff_Hz = kFreqCutoff_Hz;

    if (audioLogger.install(std::move(parameters)) == false) {
        fprintf(stderr, "Failed to install audio logger\n");
        return -1;
    }

    bool done = false;

    // Main loop
    while (!done) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            ImGui_ImplSDL2_ProcessEvent(&event);
            switch (event.type) {
                case SDL_KEYDOWN:
                    if (keyPressed == -1) {
                        //auto t1 = std::chrono::high_resolution_clock::now();
                        //printf("Event: %d\n", (int) std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count());
                        keyPressed = event.key.keysym.sym;
                        audioLogger.record(kBufferSizeTrain_s, 3);
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
        ImGui::Begin("Average Key Waveform");
        ImGui::Text("Frames in buffer: %d\n", (int) kBufferSizeTrain_frames);
        {
            std::lock_guard<std::mutex> lock(mutex);

            ImGui::PlotLines(
                "##Similarity for offset",
                similarityForOffset.data(), similarityForOffset.size(),
                0, "Similarity for offset", FLT_MAX, FLT_MAX, ImVec2(windowSizeX, 0.1f*windowSizeY));

            for (const auto & ampl : keySoundAverageAmpl) {
                struct SampleGetter {
                    static float f(void * data, int i) {
                        int fid = i/kSamplesPerFrame;
                        int sid = i - fid*kSamplesPerFrame;
                        const auto & frames = *(TKeyWaveform *)(data);
                        const auto & frame = frames[fid];
                        return frame[sid];
                    }
                };

                float (*getter)(void *, int) = SampleGetter::f;

                std::string skey(" ");
                skey[0] = ampl.first;
                ImGui::PlotLines(
                    ("##" + skey).c_str(),
                    getter, (void *)(intptr_t)(&ampl.second), ampl.second.size()*kSamplesPerFrame,
                    0, skey.c_str(), FLT_MAX, FLT_MAX, ImVec2(windowSizeX, 0.1f*windowSizeY));

                int nFrames = keySoundAverageFreq[ampl.first].size();
                for (int fid = 0; fid < (int) keySoundAverageFreq[ampl.first].size(); ++fid) {
                    const auto & freq = keySoundAverageFreq[ampl.first][fid];
                    ImGui::PlotHistogram(
                        ("##" + skey + " frame " + std::to_string(fid)).c_str(),
                        freq.data(), freq.size()/2,
                        0, ("frame " + std::to_string(fid)).c_str(), 0.0f, 50.0f, ImVec2(0.95f*windowSizeX/nFrames, 0.1f*windowSizeY));
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
