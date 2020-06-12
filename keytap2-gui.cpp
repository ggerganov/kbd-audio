/*! \file keytap2-gui.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "constants.h"
#include "common.h"
#include "common-gui.h"
#include "subbreak.h"
#include "subbreak2.h"

#include "imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl3.h"

#include <SDL.h>

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <functional>
#include <condition_variable>

static std::function<bool()> g_doInit;
static std::function<bool()> g_mainUpdate;

void mainUpdate() {
    g_mainUpdate();
}

// JS interface

extern "C" {
    int doInit() { return g_doInit(); }
}

int g_windowSizeX = 1920;
int g_windowSizeY = 1200;

struct stParameters;

using TParameters           = stParameters;

using TSampleInput          = TSampleF;
using TSample               = TSampleI16;
using TWaveform             = TWaveformI16;
using TWaveformView         = TWaveformViewI16;
using TKeyPressData         = TKeyPressDataI16;
using TKeyPressCollection   = TKeyPressCollectionI16;
using TPlaybackData         = TPlaybackDataI16;

SDL_AudioDeviceID g_deviceIdOut = 0;
TPlaybackData g_playbackData;

struct stParameters {
    int32_t playbackId              = 0;

    int32_t keyPressWidth_samples   = 256;
    int32_t offsetFromPeak_samples  = keyPressWidth_samples/2;
    int32_t alignWindow_samples     = 256;

    std::vector<int>    valuesClusters = { 50, };
    std::vector<float>  valuesPNonAlphabetic = { 0.05, 0.05, 0.05, };
    std::vector<float>  valuesWEnglishFreq = { 50.0 };

    int32_t nProcessors() const {
        return valuesClusters.size()*valuesPNonAlphabetic.size()*valuesWEnglishFreq.size();
    }

    int32_t valueForProcessorClusters(int idx) const {
      return valuesClusters[idx/(valuesPNonAlphabetic.size()*valuesWEnglishFreq.size())];
    }

    float valueForProcessorPNonAlphabetic(int idx) const {
        return valuesPNonAlphabetic[(idx/valuesWEnglishFreq.size())%valuesPNonAlphabetic.size()];
    }

    float valueForProcessorWEnglishFreq(int idx) const {
        return valuesWEnglishFreq[idx%valuesWEnglishFreq.size()];
    }

    Cipher::TParameters cipher;
};

struct stStateUI {
    bool doUpdate = false;

    struct Flags {
        bool recalculateSimilarityMap = false;
        bool resetOptimization = false;
        bool changeProcessing = false;
        bool applyClusters = false;
        bool applyWEnglishFreq = false;
        bool applyHints = false;

        void clear() { memset(this, 0, sizeof(Flags)); }
    } flags;

    bool autoHint = false;
    bool processing = true;
    bool calculatingSimilarityMap = false;

    // key presses window
    int nview = -1;
    int offset = -1;
    float amin = std::numeric_limits<TSample>::min();
    float amax = std::numeric_limits<TSample>::max();

    // window sizes
    float windowHeightKeyPesses = 0;
    float windowHeightResults = 0;
    float windowHeightSimilarity = 0;

    TParameters params;
    TWaveform waveformInput;
    TKeyPressCollection keyPresses;
    TSimilarityMap similarityMap;

    Cipher::THint suggestions;
    std::map<int, Cipher::TResult> results;
};

struct stStateCore {
    struct Flags {
        bool calculatingSimilarityMap = false;
        bool updateSimilarityMap = false;
        bool updateResult[128];

        void clear() { memset(this, 0, sizeof(Flags)); }
    } flags;

    bool processing = true;

    Cipher::TFreqMap * freqMap[3];
    TParameters params;
    TSimilarityMap similarityMap;
    TKeyPressCollection keyPresses;

    std::map<int, Cipher::Processor> processors;
};

template<typename T>
struct TripleBuffer : public T {
    bool update(bool force = false) {
        std::lock_guard<std::mutex> lock(mutex);
        hasChanged = true || force;
        buffer = *this;
        this->flags.clear();

        return true;
    }

    bool changed() {
        std::lock_guard<std::mutex> lock(mutex);
        return hasChanged;
    }

    const T & get() {
        std::lock_guard<std::mutex> lock(mutex);
        if (hasChanged) {
            out = buffer;
            hasChanged = false;
        }

        return out;
    }
private:
    bool hasChanged = false;

    T buffer;
    T out;

    mutable std::mutex mutex;
};

template<> bool TripleBuffer<stStateCore>::update(bool force) {
    std::lock_guard<std::mutex> lock(mutex);
    if (hasChanged && force == false) return false;

    hasChanged = true;

    buffer.flags = this->flags;

    if (this->flags.updateSimilarityMap) {
        buffer.similarityMap = this->similarityMap;
    }

    for (int i = 0; i < this->params.nProcessors(); ++i) {
        if (this->flags.updateResult[i]) {
            buffer.processors[i] = this->processors[i];
        }
    }

    this->flags.clear();

    return true;
}

TripleBuffer<stStateUI> stateUI;
TripleBuffer<stStateCore> stateCore;

float plotWaveform(void * data, int i) {
    TWaveformView * waveform = (TWaveformView *)data;
    return waveform->samples[i];
}

float plotWaveformInverse(void * data, int i) {
    TWaveformView * waveform = (TWaveformView *)data;
    return -waveform->samples[i];
}

bool renderKeyPresses(stStateUI & stateUI, const char * fnameInput, const TWaveform & waveform, TKeyPressCollection & keyPresses) {
    auto & params = stateUI.params;

    int & nview = stateUI.nview;
    int & offset = stateUI.offset;
    float & amin = stateUI.amin;
    float & amax = stateUI.amax;

    if (nview < 0) nview = waveform.size();
    if (offset < 0) offset = (waveform.size() - nview)/2;

    ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Once);
    ImGui::SetNextWindowSize(ImVec2(g_windowSizeX, 340.0f), ImGuiCond_Always);
    if (ImGui::Begin("Key Presses", nullptr, ImGuiWindowFlags_NoMove)) {
        int viewMin = 512;
        int viewMax = waveform.size();

        bool ignoreDelete = false;

        static float dragOffset = 0.0f;
        static float scrollSize = 18.0f;

        static auto nviewPrev = nview + 1;

        static TWaveform waveformLowRes = waveform;
        static TWaveform waveformThreshold = waveform;
        static TWaveform waveformMax = waveform;

        auto wview = getView(waveformLowRes, offset, nview);
        auto tview = getView(waveformThreshold, offset, nview);
        auto mview = getView(waveformMax, offset, nview);
        auto wsize = ImVec2(ImGui::GetContentRegionAvailWidth(), 250.0);

        auto mpos = ImGui::GetIO().MousePos;
        auto savePos = ImGui::GetCursorScreenPos();
        auto drawList = ImGui::GetWindowDrawList();
        ImGui::PushStyleColor(ImGuiCol_FrameBg, { 0.3f, 0.3f, 0.3f, 0.3f });
        ImGui::PushStyleColor(ImGuiCol_PlotHistogram, { 1.0f, 1.0f, 1.0f, 1.0f });
        ImGui::PlotHistogram("##Waveform", plotWaveform, &wview, nview, 0, "Waveform", amin, amax, wsize);
        ImGui::PopStyleColor(2);
        ImGui::SetCursorScreenPos(savePos);
        ImGui::PushStyleColor(ImGuiCol_FrameBg, { 0.1f, 0.1f, 0.1f, 0.0f });
        ImGui::PushStyleColor(ImGuiCol_PlotHistogram, { 1.0f, 1.0f, 1.0f, 1.0f });
        ImGui::PlotHistogram("##Waveform", plotWaveformInverse, &wview, nview, 0, "Waveform", amin, amax, wsize);
        ImGui::PopStyleColor(2);
        ImGui::SetCursorScreenPos(savePos);
        ImGui::PushStyleColor(ImGuiCol_FrameBg, { 0.0f, 0.0f, 0.0f, 0.0f });
        ImGui::PushStyleColor(ImGuiCol_PlotLines, { 1.0f, 1.0f, 0.0f, 0.5f });
        ImGui::PlotLines("##WaveformThreshold", plotWaveform, &tview, nview, 0, "", amin, amax, wsize);
        ImGui::PopStyleColor(2);
        ImGui::SetCursorScreenPos(savePos);
        ImGui::PushStyleColor(ImGuiCol_FrameBg, { 0.0f, 0.0f, 0.0f, 0.0f });
        ImGui::PushStyleColor(ImGuiCol_PlotLines, { 1.0f, 1.0f, 0.0f, 0.5f });
        ImGui::PlotLines("##WaveformThreshold", plotWaveformInverse, &tview, nview, 0, "", amin, amax, wsize);
        ImGui::PopStyleColor(2);
        ImGui::SetCursorScreenPos(savePos);
        ImGui::PushStyleColor(ImGuiCol_FrameBg, { 0.0f, 0.0f, 0.0f, 0.0f });
        ImGui::PushStyleColor(ImGuiCol_PlotLines, { 0.0f, 1.0f, 0.0f, 0.5f });
        ImGui::PlotLines("##WaveformMax", plotWaveform, &mview, nview, 0, "", amin, amax, wsize);
        ImGui::PopStyleColor(2);
        ImGui::SetCursorScreenPos(savePos);
        ImGui::PushStyleColor(ImGuiCol_FrameBg, { 0.0f, 0.0f, 0.0f, 0.0f });
        ImGui::PushStyleColor(ImGuiCol_PlotLines, { 0.0f, 1.0f, 0.0f, 0.5f });
        ImGui::PlotLines("##WaveformMax", plotWaveformInverse, &mview, nview, 0, "", amin, amax, wsize);
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

            if (ImGui::IsMouseReleased(0) && ImGui::GetIO().KeyCtrl) {
                int i = 0;
                int64_t pos = offset + nview*(mpos.x - savePos.x)/wsize.x;
                for (i = 0; i < (int) keyPresses.size(); ++i) {
                    int64_t p0 = keyPresses[i].pos + params.offsetFromPeak_samples - params.keyPressWidth_samples;
                    int64_t p1 = keyPresses[i].pos + params.offsetFromPeak_samples + params.keyPressWidth_samples;
                    int64_t pmin = std::min(std::min(keyPresses[i].pos, p0), p1);
                    int64_t pmax = std::max(std::max(keyPresses[i].pos, p0), p1);
                    if (pmin > pos) {
                        ignoreDelete = true;
                        TKeyPressData entry;
                        entry.pos = pos;
                        entry.waveform = getView(waveform, 0);
                        keyPresses.insert(keyPresses.begin() + i, entry);
                        break;
                    }
                    if (pos >= pmin && pos <= pmax) break;
                }
                if (i == (int) keyPresses.size() && ignoreDelete == false) {
                    ignoreDelete = true;
                    TKeyPressData entry;
                    entry.pos = pos;
                    entry.waveform = getView(waveform, 0);
                    keyPresses.push_back(entry);
                }
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
        for (int i = 0; i < (int) keyPresses.size(); ++i) {
            if (keyPresses[i].pos + params.offsetFromPeak_samples + params.keyPressWidth_samples < offset) continue;
            if (keyPresses[i].pos + params.offsetFromPeak_samples - params.keyPressWidth_samples >= offset + nview) break;

            {
                float x0 = ((float)(keyPresses[i].pos - offset))/nview;

                ImVec2 p0 = { savePos.x + x0*wsize.x, savePos.y };
                ImVec2 p1 = { savePos.x + x0*wsize.x, savePos.y + wsize.y };

                drawList->AddLine(p0, p1, ImGui::ColorConvertFloat4ToU32({ 1.0f, 0.0f, 0.0f, 0.8f }), 1.0f);
            }

            {
                float x0 = ((float)(keyPresses[i].pos - offset))/nview;
                float x1 = ((float)(keyPresses[i].pos + params.offsetFromPeak_samples - params.keyPressWidth_samples - offset))/nview;
                float x2 = ((float)(keyPresses[i].pos + params.offsetFromPeak_samples + params.keyPressWidth_samples - offset))/nview;

                ImVec2 p0 = { savePos.x + x0*wsize.x, savePos.y };
                ImVec2 p1 = { savePos.x + x1*wsize.x, savePos.y };
                ImVec2 p2 = { savePos.x + x2*wsize.x, savePos.y + wsize.y };

                float xmin = std::min(std::min(p0.x, p1.x), p2.x);
                float xmax = std::max(std::max(p0.x, p1.x), p2.x);

                bool isHovered = (mpos.x > xmin && mpos.x < xmax && mpos.y > p1.y && mpos.y < p2.y);

                float col = isHovered ? 0.7f : 0.3f;
                drawList->AddRectFilled(p1, p2, ImGui::ColorConvertFloat4ToU32({ 1.0f, 0.0f, 0.0f, col }));

                if (isHovered) {
                    if (ImGui::IsMouseReleased(0) && ImGui::GetIO().KeyCtrl && ignoreDelete == false) {
                        keyPresses.erase(keyPresses.begin() + i);
                        --i;
                    }
                }

                if (nview < 64.0*wsize.x) {
                    ImGui::SetCursorScreenPos({ savePos.x + 0.5f*((x1 + x2)*wsize.x - ImGui::CalcTextSize(std::to_string(i).c_str()).x), savePos.y + wsize.y - ImGui::GetTextLineHeightWithSpacing() });
                    ImGui::Text("%d", i);
                }

            }

            {
                float x0 = ((float)(g_playbackData.offset + g_playbackData.idx - offset))/nview;

                ImVec2 p0 = {savePos.x + x0*wsize.x, savePos.y};
                ImVec2 p1 = {savePos.x + x0*wsize.x, savePos.y + wsize.y};
                drawList->AddLine(p0, p1, ImGui::ColorConvertFloat4ToU32({ 1.0f, 1.0f, 0.0f, 0.3f }));
            }
        }
        for (int i = 0; i < (int) keyPresses.size(); ++i) {
            {
                float x0 = ((float)(keyPresses[i].pos))/viewMax;

                ImVec2 p0 = { savePos.x + x0*wsize.x, savePos.y + wsize.y };
                ImVec2 p1 = { savePos.x + x0*wsize.x, savePos.y + wsize.y + scrollSize };

                drawList->AddLine(p0, p1, ImGui::ColorConvertFloat4ToU32({ 1.0f, 0.0f, 0.0f, 0.3f }));
            }
        }

        ImGui::SetCursorScreenPos(savePos2);

        //auto io = ImGui::GetIO();
        //ImGui::Text("Keys pressed:");   for (int i = 0; i < IM_ARRAYSIZE(io.KeysDown); i++) if (ImGui::IsKeyPressed(i))             { ImGui::SameLine(); ImGui::Text("%d", i); }

        static bool recalculate = true;
        static bool playHalfSpeed = false;
        static int historySize = 3*1024;
        static float thresholdBackground = 10.0;
        ImGui::PushItemWidth(100.0);

        ImGui::Checkbox("x0.5", &playHalfSpeed);
        ImGui::SameLine();
        if (g_playbackData.playing) {
            if (ImGui::Button("Stop") || (ImGui::IsKeyPressed(44) && ImGui::GetIO().KeyCtrl)) { // Ctrl + space
                g_playbackData.playing = false;
                g_playbackData.idx = g_playbackData.waveform.n - TPlaybackData::kSamples;
            }
        } else {
            if (ImGui::Button("Play") || (ImGui::IsKeyPressed(44) && ImGui::GetIO().KeyCtrl)) { // Ctrl + space
                g_playbackData.playing = true;
                g_playbackData.slowDown = playHalfSpeed ? 2 : 1;
                g_playbackData.idx = 0;
                g_playbackData.offset = offset;
                g_playbackData.waveform = getView(waveform, offset, std::min((int) (10*kSampleRate), nview));
                SDL_PauseAudioDevice(g_deviceIdOut, 0);
            }
        }

        if (g_playbackData.idx >= g_playbackData.waveform.n) {
            g_playbackData.playing = false;
            SDL_ClearQueuedAudio(g_deviceIdOut);
            SDL_PauseAudioDevice(g_deviceIdOut, 1);
        }

        ImGui::SameLine();
        ImGui::SliderFloat("Threshold background", &thresholdBackground, 0.1f, 50.0f) && (recalculate = true);
        ImGui::SameLine();
        ImGui::SliderInt("History Size", &historySize, 512, 1024*16) && (recalculate = true);
        ImGui::SameLine();
        if (ImGui::Button("Recalculate") || recalculate) {
            findKeyPresses(getView(waveform, 0), keyPresses, waveformThreshold, waveformMax, thresholdBackground, historySize, true);
            recalculate = false;
        }

        static std::string filename = std::string(fnameInput) + ".keys";
        ImGui::SameLine();
        ImGui::Text("%s", filename.c_str());

        ImGui::SameLine();
        if (ImGui::Button("Save")) {
            saveKeyPresses(filename.c_str(), keyPresses);
        }

        ImGui::SameLine();
        if (ImGui::Button("Load")) {
            loadKeyPresses(filename.c_str(), getView(waveform, 0), keyPresses);
        }

        ImGui::SameLine();
        ImGui::DragInt("Key width", &params.keyPressWidth_samples, 8, 0, kSampleRate/10);
        ImGui::SameLine();
        ImGui::DragInt("Peak offset", &params.offsetFromPeak_samples, 8, -kSampleRate/10, kSampleRate/10);
        ImGui::SameLine();
        ImGui::DragInt("Align window", &params.alignWindow_samples, 8, 0, kSampleRate/10);

        ImGui::PopItemWidth();

        if (nview != nviewPrev) {
            generateLowResWaveform(waveform, waveformLowRes, std::max(1.0f, nview/wsize.x));
            nviewPrev = nview;
        }

        stateUI.windowHeightKeyPesses = ImGui::GetWindowHeight();
    } else {
        stateUI.windowHeightKeyPesses = ImGui::GetTextLineHeightWithSpacing();
    }
    ImGui::End();

    return false;
}

bool renderResults(stStateUI & stateUI) {
    ImGui::SetNextWindowPos(ImVec2(0, stateUI.windowHeightKeyPesses), ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(1.0f*g_windowSizeX, (12 + stateUI.results.size())*ImGui::GetTextLineHeightWithSpacing()), ImGuiCond_Always);
    if (ImGui::Begin("Results", nullptr, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_NoMove)) {
        auto drawList = ImGui::GetWindowDrawList();

        if (ImGui::Button("Calculate Key Similarity")) {
            stateUI.similarityMap.clear();
            stateUI.flags.recalculateSimilarityMap = true;
            stateUI.doUpdate = true;
        }

        if (stateUI.similarityMap.size() > 0 && stateUI.calculatingSimilarityMap == false) {
            ImGui::SameLine();
            if (stateUI.processing) {
                if (ImGui::Button("Pause") || (ImGui::IsKeyPressed(44) && !ImGui::GetIO().KeyCtrl)) { // space
                    stateUI.processing = false;
                    stateUI.flags.changeProcessing = true;
                    stateUI.doUpdate = true;
                }
            } else {
                if (ImGui::Button("Resume") || (ImGui::IsKeyPressed(44) && !ImGui::GetIO().KeyCtrl)) { // space
                    stateUI.processing = true;
                    stateUI.flags.changeProcessing = true;
                    stateUI.doUpdate = true;
                }
            }

            ImGui::SameLine();
            if (ImGui::Button("Reset")) { // c
                stateUI.flags.resetOptimization = true;
                stateUI.doUpdate = true;
            }

            ImGui::SameLine();
            if (ImGui::Button("Apply Suggestions")) {
                int n = stateUI.suggestions.size();
                for (int i = 0; i < n; ++i) {
                    if (stateUI.suggestions[i] < 0) {
                        stateUI.keyPresses[i].bind = -1;
                        continue;
                    }
                    if (stateUI.suggestions[i] > 0 && stateUI.suggestions[i] <= 26) {
                        stateUI.keyPresses[i].bind = stateUI.suggestions[i] - 1;
                    } else {
                        stateUI.keyPresses[i].bind = 26;
                    }
                }

                stateUI.flags.applyHints = true;
                stateUI.doUpdate = true;
            }

            ImGui::SameLine();
            if (ImGui::Button("Clear Hints")) {
                for (auto & keyPress : stateUI.keyPresses) {
                    keyPress.bind = -1;
                }
                stateUI.flags.applyHints = true;
                stateUI.doUpdate = true;
            }
        } else {
            if (stateUI.calculatingSimilarityMap) {
                ImGui::SameLine();
                ImGui::TextColored({0.3f, 1.0f, 0.4f, 1.0f}, "Calculating similarity map ...");
            }
        }

        ImGui::PushItemWidth(200.0);

        if (ImGui::Checkbox("Auto", &stateUI.autoHint)) {
        }
        if (ImGui::IsItemHovered()) {
            ImGui::BeginTooltip();
            ImGui::Text("If checked, hints will be updated automatically based on the suggestions");
            ImGui::EndTooltip();
        }

        if (stateUI.autoHint && stateUI.results.size() > 0) {
            int nPerAutoHint = 20;
            auto id0 = stateUI.results.begin()->second.id%stateUI.params.cipher.nMHIters;
            auto id1 = (stateUI.results.begin()->second.id/stateUI.params.cipher.nMHIters)%nPerAutoHint;
            ImGui::SameLine();
            ImGui::ProgressBar(float(id1)/nPerAutoHint, { 100.0, ImGui::GetTextLineHeightWithSpacing() });
            if (id0 == 0 && id1 == 0) {
                int n = stateUI.suggestions.size();
                for (int i = 0; i < n; ++i) {
                    if (frand() > 0.1) {
                        stateUI.keyPresses[i].bind = -1;
                        continue;
                    }
                    if (stateUI.suggestions[i] < 0) {
                        stateUI.keyPresses[i].bind = -1;
                        continue;
                    }
                    //if (stateUI.suggestions[i] == 0 || stateUI.suggestions[i] == 5) {
                    //    stateUI.keyPresses[i].bind = frand() > 0.5 ? 4 : 26;
                    //    continue;
                    //}
                    if (stateUI.suggestions[i] > 0 && stateUI.suggestions[i] <= 26) {
                        stateUI.keyPresses[i].bind = stateUI.suggestions[i] - 1;
                    } else {
                        stateUI.keyPresses[i].bind = 26;
                    }
                }
                stateUI.flags.applyHints = true;
                stateUI.flags.applyWEnglishFreq = true;
                stateUI.doUpdate = true;
            }
        }

        ImGui::SameLine();
        if (ImGui::SliderInt("Clusters", &stateUI.params.valuesClusters[0], 30, 150.0)) {
            stateUI.flags.applyClusters = true;
            stateUI.doUpdate = true;
        }

        //ImGui::SameLine();
        //if (ImGui::SliderFloat("English weight", &stateUI.params.valuesWEnglishFreq[0], 1.0, 100.0, "%.6f", 2.0)) {
        //    stateUI.flags.applyWEnglishFreq = true;
        //    stateUI.doUpdate = true;
        //}

        ImGui::PopItemWidth();

        ImGui::Separator();
        ImGui::Text(" Suggestions");
        {
            int n = stateUI.suggestions.size();
            if (n > 0) {
                ImGui::SameLine();
                auto charSize = ImGui::CalcTextSize("a");
                auto curPos = ImGui::GetCursorScreenPos();
                ImGui::InvisibleButton("", ImVec2{charSize.x*n, 1.0f});
                ImGui::SetCursorScreenPos(curPos);
                for (int i = 0; i < n; ++i) {
                    char c = ' ';
                    if (stateUI.suggestions[i] < 0) {
                        ImGui::Text(" ");
                    } else {
                        c = '_';
                        if (stateUI.suggestions[i] > 0 && stateUI.suggestions[i] <= 26) {
                            c = 'a'+stateUI.suggestions[i] - 1;
                        }
                        ImGui::TextColored({ 1.0f, 1.0f, 0.0f, 1.0f }, "%c", c);
                    }
                    if (ImGui::IsItemHovered()) {
                        ImGui::BeginTooltip();
                        ImGui::Text("Suggestion for key press %d: '%c'", i, c);
                        ImGui::EndTooltip();
                    }
                    if (i < n - 1) {
                        curPos.x += charSize.x;
                        ImGui::SetCursorScreenPos(curPos);
                    }
                }
            }
        }

        ImGui::Text("       Hints");
        {
            int n = stateUI.keyPresses.size();
            if (n > 0) {
                ImGui::SameLine();
                auto charSize = ImGui::CalcTextSize("a");
                auto curPos = ImGui::GetCursorScreenPos();
                ImGui::InvisibleButton("", ImVec2{charSize.x*n, 1.0f});
                ImGui::SetCursorScreenPos(curPos);
                for (int i = 0; i < n; ++i) {
                    if (stateUI.keyPresses[i].bind < 0) {
                        ImGui::Text(".");
                    } else {
                        char c = '_';
                        if (stateUI.keyPresses[i].bind >= 0 && stateUI.keyPresses[i].bind < 26) {
                            c = 'a'+stateUI.keyPresses[i].bind;
                        }
                        ImGui::TextColored({ 0.0f, 1.0f, 0.0f, 1.0f }, "%c", c);
                    }
                    if (ImGui::IsItemHovered()) {
                        auto p0 = curPos;
                        auto p1 = p0;
                        p1.x += charSize.x;
                        p1.y += charSize.y;
                        drawList->AddRect(p0, p1, ImGui::ColorConvertFloat4ToU32({0.0f, 1.0f, 0.0f, 1.0f}));
                        drawList->AddRectFilled(p0, p1, ImGui::ColorConvertFloat4ToU32({0.0f, 1.0f, 0.0f, 0.6f}));

                        if (ImGui::IsMouseDown(0)) {
                            stateUI.nview = 16*1024;
                            stateUI.offset = stateUI.keyPresses[i].pos - stateUI.nview/2;
                        }

                        if (ImGui::IsMouseDown(1)) {
                            stateUI.keyPresses[i].bind = -1;
                            stateUI.flags.applyHints = true;
                            stateUI.doUpdate = true;
                        }

                        for (int k = 0; k < 26; ++k) {
                            if (ImGui::IsKeyPressed(4 + k)) {
                                stateUI.keyPresses[i].bind = k;
                            }
                        }

                        ImGui::BeginTooltip();
                        ImGui::Text("Key press: %d", i);
                        ImGui::Text(" - Left click to focus");
                        ImGui::Text(" - Right click to remove hint");
                        ImGui::EndTooltip();
                    }
                    if (i < n - 1) {
                        curPos.x += charSize.x;
                        ImGui::SetCursorScreenPos(curPos);
                    }
                }
            }
        }

        ImGui::Text("%2s. %8s", "#", "Prob");
        ImGui::Separator();

        for (const auto & item : stateUI.results) {
            const auto & id = item.first;
            const auto & result = item.second;
            ImGui::PushID(id);
            ImGui::Text("%2d. %8.3f", id, result.p);
            if (ImGui::IsItemHovered()) {
                ImGui::BeginTooltip();
                ImGui::Text("Processor   %d", item.first);
                ImGui::Text("Cluster-to-letter map:");
                {
                    int i = 0;
                    for (auto& pair : item.second.clMap) {
                        ImGui::Text("%3d: %3d  ", pair.first, pair.second);
                        if (++i % 10 != 0) ImGui::SameLine();
                    }
                }
                ImGui::EndTooltip();
            }
            if (result.clMap.empty() == false) {
                ImGui::SameLine();
                int n = result.clusters.size();
                auto charSize = ImGui::CalcTextSize("a");
                auto curPos = ImGui::GetCursorScreenPos();
                ImGui::InvisibleButton("", ImVec2{charSize.x*n, 1.0f});
                ImGui::SetCursorScreenPos(curPos);
                for (int i = 0; i < n; ++i) {
                    auto cluster = result.clusters[i];
                    if (result.clMap.find(cluster) == result.clMap.end()) {
                        ImGui::Text("%c", '?');
                    } else {
                        auto let = result.clMap.at(result.clusters[i]);

                        char c = '.';
                        if (let > 0 && let <= 26) {
                            c = 'a'+let - 1;
                        }

                        if (stateUI.keyPresses[i].bind == let - 1) {
                            ImGui::TextColored({ 0.0f, 1.0f, 0.0f, 1.0f }, "%c", c);
                        } else if (stateUI.suggestions[i] == let) {
                            ImGui::TextColored({ 1.0f, 1.0f, 0.0f, 1.0f }, "%c", c);
                        } else {
                            ImGui::Text("%c", c);
                        }

                        if (ImGui::IsItemHovered()) {
                            auto p0 = curPos;
                            auto p1 = p0;
                            p1.x += charSize.x;
                            p1.y += charSize.y;
                            drawList->AddRect(p0, p1, ImGui::ColorConvertFloat4ToU32({0.0f, 1.0f, 0.0f, 1.0f}));
                            drawList->AddRectFilled(p0, p1, ImGui::ColorConvertFloat4ToU32({0.0f, 1.0f, 0.0f, 0.6f}));
                            if (ImGui::IsMouseDown(0)) {
                                if (let > 0 && let <= 26) {
                                    stateUI.keyPresses[i].bind = let - 1;
                                } else {
                                    stateUI.keyPresses[i].bind = 26;
                                }
                                stateUI.flags.applyHints = true;
                                stateUI.doUpdate = true;
                            }
                            ImGui::BeginTooltip();
                            ImGui::Text("Key presss: %d, letter - %d", i, let);
                            ImGui::Text(" - Left click to set as hint");
                            ImGui::EndTooltip();
                        }
                    }

                    if (i < n - 1) {
                        curPos.x += charSize.x;
                        ImGui::SetCursorScreenPos(curPos);
                    }
                }
            }
            ImGui::PopID();
        }

        stateUI.windowHeightResults = ImGui::GetWindowHeight();
    } else {
        stateUI.windowHeightResults = ImGui::GetTextLineHeightWithSpacing();
    }

    ImGui::End();
    return true;
}

bool renderSimilarity(const TKeyPressCollection & keyPresses, const TSimilarityMap & similarityMap) {
    int offsetY = stateUI.windowHeightKeyPesses + stateUI.windowHeightResults;
    ImGui::SetNextWindowPos(ImVec2(0, offsetY), ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(1.0f*g_windowSizeX, g_windowSizeY - offsetY), ImGuiCond_Always);
    if (ImGui::Begin("Similarity", nullptr, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_NoMove)) {
        auto wsize = ImGui::GetContentRegionAvail();

        static float bsize = 4.0f;
        static float threshold = 0.0f;
        ImGui::PushItemWidth(100.0);

        int n = similarityMap.size();
        ImGui::SameLine();
        ImGui::SliderFloat("Size", &bsize, 1.5f, 24.0f);
        ImGui::SameLine();
        if (ImGui::Button("Fit")) {
            bsize = (std::min(wsize.x, wsize.y) - 24.0)/n;
        }
        ImGui::SameLine();
        ImGui::SliderFloat("Threshold", &threshold, 0.0f, 1.0f);
        ImGui::PopItemWidth();

        ImGui::BeginChild("Canvas", { 0, 0 }, 1, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_AlwaysHorizontalScrollbar | ImGuiWindowFlags_AlwaysVerticalScrollbar);
        auto savePos = ImGui::GetCursorScreenPos();
        auto drawList = ImGui::GetWindowDrawList();

        int hoveredId = -1;
        if (n > 0) {
            ImGui::InvisibleButton("SimilarityMapIB", { n*bsize, n*bsize });
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    float col = similarityMap[i][j].cc;
                    ImVec2 p0 = {savePos.x + j*bsize, savePos.y + i*bsize};
                    ImVec2 p1 = {savePos.x + (j + 1)*bsize - 1.0f, savePos.y + (i + 1)*bsize - 1.0f};
                    if (similarityMap[i][j].cc > threshold) {
                        drawList->AddRectFilled(p0, p1, ImGui::ColorConvertFloat4ToU32(ImVec4{1.0f - col, col, 0.0f, col}));
                    }
                    if (ImGui::IsMouseHoveringRect(p0, {p1.x + 1.0f, p1.y + 1.0f})) {
                        if (i == j) hoveredId = i;
                        if (ImGui::IsMouseDown(0) == false) {
                            ImGui::BeginTooltip();
                            ImGui::Text("[%3d, %3d]\n", keyPresses[i].cid, keyPresses[j].cid);
                            ImGui::Text("[%3d, %3d] = %5.4g\n", i, j, similarityMap[i][j].cc);
                            for (int k = 0; k < n; ++k) {
                                if (similarityMap[i][k].cc > 0.5) ImGui::Text("Offset [%3d, %3d] = %d\n", i, k, (int) similarityMap[i][k].offset);
                            }
                            ImGui::Separator();
                            for (int k = 0; k < n; ++k) {
                                if (similarityMap[k][i].cc > 0.5) ImGui::Text("Offset [%3d, %3d] = %d\n", k, i, (int) similarityMap[k][i].offset);
                            }
                            ImGui::EndTooltip();
                        }
                    }
                }
            }

            if (hoveredId >= 0) {
                for (int i = 0; i < n; ++i) {
                    if (keyPresses[i].cid == keyPresses[hoveredId].cid) {
                        ImVec2 p0 = {savePos.x + i*bsize, savePos.y + i*bsize};
                        ImVec2 p1 = {savePos.x + (i + 1)*bsize - 1.0f, savePos.y + (i + 1)*bsize - 1.0f};
                        drawList->AddRectFilled(p0, p1, ImGui::ColorConvertFloat4ToU32({1.0f, 0.0f, 0.0f, 1.0f}));
                    }
                }
            }
        }

        ImGui::EndChild();
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

    playbackSpec.freq = kSampleRate;
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
    srand(time(0));

    printf("Usage: %s recrod.kbd n-gram-dir [letter.mask] [-pN]\n", argv[0]);
    printf("    -pN - select playback device N\n");
    if (argc < 3) {
        return -1;
    }

    auto argm = parseCmdArguments(argc, argv);
    int playbackId = argm["p"].empty() ? 0 : std::stoi(argm["p"]);

    stateUI.params.playbackId = playbackId;

    if (SDL_Init(SDL_INIT_VIDEO|SDL_INIT_TIMER) != 0) {
        printf("Error: %s\n", SDL_GetError());
        return -2;
    }

    g_doInit = [&]() {
        if (prepareAudioOut(stateUI.params) == false) {
            printf("Error: failed to initialize audio playback\n");
            return false;
        }

        return true;
    };

    printf("[+] Loading recording from '%s'\n", argv[1]);
    if (readFromFile<TSampleF>(argv[1], stateUI.waveformInput) == false) {
        printf("Specified file '%s' does not exist\n", argv[1]);
        return -4;
    }

    Cipher::TFreqMap freqMap3;
    Cipher::TFreqMap freqMap4;
    Cipher::TFreqMap freqMap5;

    if (Cipher::loadFreqMap((std::string(argv[2]) + "/./english_trigrams.txt").c_str(), freqMap3) == false) {
        return -5;
    }
    if (Cipher::loadFreqMap((std::string(argv[2]) + "/./english_quadgrams.txt").c_str(), freqMap4) == false) {
        return -5;
    }
    if (Cipher::loadFreqMap((std::string(argv[2]) + "/./english_quintgrams.txt").c_str(), freqMap5) == false) {
        return -5;
    }
    stateCore.freqMap[0] = &freqMap3;
    stateCore.freqMap[1] = &freqMap4;
    stateCore.freqMap[2] = &freqMap5;

    Gui::Objects guiObjects;
    if (Gui::init("Keytap2", g_windowSizeX, g_windowSizeY, guiObjects) == false) {
        return -6;
    }

    printf("[+] Loaded recording: of %d samples (sample size = %d bytes)\n", (int) stateUI.waveformInput.size(), (int) sizeof(TSample));
    printf("    Size in memory:          %g MB\n", (float)(sizeof(TSample)*stateUI.waveformInput.size())/1024/1024);
    printf("    Sample size:             %d\n", (int) sizeof(TSample));
    printf("    Total number of samples: %d\n", (int) stateUI.waveformInput.size());
    printf("    Recording length:        %g seconds\n", (float)(stateUI.waveformInput.size())/kSampleRate);

    bool finishApp = false;
    g_mainUpdate = [&]() {
        if (finishApp) return false;

        if (stateCore.changed()) {
            auto stateCoreNew = stateCore.get();

            stateUI.calculatingSimilarityMap = stateCoreNew.flags.calculatingSimilarityMap;

            if (stateCoreNew.flags.updateSimilarityMap) {
                stateUI.similarityMap = stateCoreNew.similarityMap;
            }

            bool recalcSuggestions = false;
            for (int i = 0; i < stateCoreNew.params.nProcessors(); ++i) {
                if (stateCoreNew.flags.updateResult[i]) {
                    recalcSuggestions = true;
                    stateUI.results[i] = stateCoreNew.processors[i].getResult();
                }
            }

            if (recalcSuggestions) {
                int n = stateUI.keyPresses.size();
                stateUI.suggestions.resize(n);

                for (int i = 0; i < n; ++i) {
                    int c = 0;
                    std::map<int, int> cnt;
                    for (const auto & result : stateUI.results) {
                        if (result.second.clMap.find(result.second.clusters[i]) == result.second.clMap.end()) continue;
                        if (i <= (int) result.second.clusters.size()) {
                            ++cnt[result.second.clMap.at(result.second.clusters[i])];
                            ++c;
                        }
                    }
                    bool hasSuggestion = false;
                    for (auto & s : cnt) {
                        if (s.second > c/2) {
                            stateUI.suggestions[i] = s.first;
                            hasSuggestion = true;
                        }
                    }
                    if (hasSuggestion == false) {
                        stateUI.suggestions[i] = -1;
                    }
                }
            }
        }

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
                    if (event.window.event == SDL_WINDOWEVENT_CLOSE && event.window.windowID == SDL_GetWindowID(guiObjects.window)) {
                        finishApp = true;
                    }
                    break;
            };
        }

        SDL_GetWindowSize(guiObjects.window, &g_windowSizeX, &g_windowSizeY);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplSDL2_NewFrame(guiObjects.window);
        ImGui::NewFrame();

        for (int i = 0; i < IM_ARRAYSIZE(ImGui::GetIO().KeysDown); i++) if (ImGui::IsKeyPressed(i)) {
            //printf("Key pressed: %d\n", i);
        }

        renderKeyPresses(stateUI, argv[1], stateUI.waveformInput, stateUI.keyPresses);
        renderResults(stateUI);
        renderSimilarity(stateUI.keyPresses, stateUI.similarityMap);

        Gui::render(guiObjects);

        if (stateUI.doUpdate) {
            stateUI.update();
            stateUI.doUpdate = false;
        }

        return true;
    };

    std::thread worker([&]() {
        while (finishApp == false) {
            if (stateUI.changed()) {
                auto stateUINew = stateUI.get();

                if (stateUINew.flags.recalculateSimilarityMap || stateUINew.flags.resetOptimization) {
                    printf("[+] Recalculating similarity map ...\n");

                    stateCore.params = stateUINew.params;
                    stateCore.keyPresses = stateUINew.keyPresses;

                    if (stateUINew.flags.recalculateSimilarityMap) {
                        stateCore.flags.calculatingSimilarityMap = true;
                        stateCore.update(true);

                        calculateSimilartyMap(
                                stateUINew.params.keyPressWidth_samples,
                                stateUINew.params.alignWindow_samples,
                                stateUINew.params.offsetFromPeak_samples,
                                stateCore.keyPresses,
                                stateCore.similarityMap);

                        if (adjustKeyPresses(stateCore.keyPresses, stateCore.similarityMap)) {
                            calculateSimilartyMap(
                                    stateUINew.params.keyPressWidth_samples,
                                    stateUINew.params.alignWindow_samples,
                                    stateUINew.params.offsetFromPeak_samples,
                                    stateCore.keyPresses,
                                    stateCore.similarityMap);
                        }

                        printf("[+] Similarity map recalculated\n");

                        stateCore.flags.calculatingSimilarityMap = false;
                        stateCore.flags.updateSimilarityMap = true;
                        stateCore.update();
                    }

                    int n = stateCore.params.nProcessors();
                    for (int i = 0; i < n; ++i) {
                        int nClusters = stateCore.params.valueForProcessorClusters(i);
                        double p = stateCore.params.valueForProcessorPNonAlphabetic(i);
                        double w = stateCore.params.valueForProcessorWEnglishFreq(i);

                        Cipher::TParameters params;
                        params.maxClusters = nClusters;
                        stateCore.processors[i].init(
                                params,
                                *stateCore.freqMap[i%3],
                                stateCore.similarityMap);

                        stateCore.processors[i].setPNonAlphabetic(std::log(p));
                        stateCore.processors[i].setWEnglishFreq(w);

                        printf("[+] Processor %d initialized: cluster = %d, p = %g, w = %g\n", i, nClusters, p, w);
                    }
                }

                if (stateUINew.flags.changeProcessing) {
                    stateCore.processing = stateUINew.processing;
                }

                if (stateUINew.flags.applyClusters) {
                    stateCore.params = stateUINew.params;

                    int n = stateCore.params.nProcessors();
                    for (int i = 0; i < n; ++i) {
                        int nClusters = stateCore.params.valueForProcessorClusters(i);
                        double p = stateCore.params.valueForProcessorPNonAlphabetic(i);
                        double w = stateCore.params.valueForProcessorWEnglishFreq(i);

                        Cipher::TParameters params;
                        params.maxClusters = nClusters;
                        stateCore.processors[i].init(
                                params,
                                *stateCore.freqMap[i%3],
                                stateCore.similarityMap);

                        stateCore.processors[i].setPNonAlphabetic(std::log(p));
                        stateCore.processors[i].setWEnglishFreq(w);
                    }
                }

                if (stateUINew.flags.applyWEnglishFreq) {
                    stateCore.params = stateUINew.params;

                    int n = stateCore.params.nProcessors();
                    for (int i = 0; i < n; ++i) {
                        double w = stateCore.params.valueForProcessorWEnglishFreq(i);
                        stateCore.processors[i].setWEnglishFreq(w);
                    }
                }

                if (stateUINew.flags.applyHints) {
                    stateCore.keyPresses = stateUINew.keyPresses;
                }
            }

            if (stateCore.processing && stateCore.keyPresses.size() > 0) {
                {
                    int n = stateCore.keyPresses.size();
                    stateCore.params.cipher.hint.clear();
                    stateCore.params.cipher.hint.resize(n, -1);
                    for (int i = 0; i < n; ++i) {
                        if (stateCore.keyPresses[i].bind < 0) continue;
                        stateCore.params.cipher.hint[i] = stateCore.keyPresses[i].bind + 1;
                    }
                }

                int nFinished = 0;
                int nWorkers = std::min(stateCore.params.nProcessors(), (int) std::thread::hardware_concurrency()/2);

                std::mutex mutex;
                std::condition_variable cv;
                std::vector<std::thread> workers(nWorkers);
                for (int iw = 0; iw < (int) workers.size(); ++iw) {
                    auto & worker = workers[iw];
                    worker = std::thread([&](int ith) {
                        int n = stateCore.params.nProcessors();
                        for (int i = ith; i < n; i += nWorkers) {
                            stateCore.processors[i].setHint(stateCore.params.cipher.hint);
                            stateCore.processors[i].compute();

                            stateCore.flags.updateResult[i] = true;
                            stateCore.update();
                        }
                        {
                            std::lock_guard<std::mutex> lock(mutex);
                            ++nFinished;
                            cv.notify_one();
                        }
                    }, iw);
                    worker.detach();
                }

                std::unique_lock<std::mutex> lock(mutex);
                cv.wait(lock, [&]() { return nFinished == nWorkers; });
            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
    });

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

    worker.join();

    printf("[+] Terminated\n");

    Gui::free(guiObjects);

    return 0;
}
