#ifdef __EMSCRIPTEN__
#include "emscripten.h"
#else
#define EMSCRIPTEN_KEEPALIVE
#endif

#include "build-vars.h"
#include "constants.h"
#include "common.h"
#include "common-gui.h"
#include "subbreak3.h"
#include "audio-logger.h"

#include "imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl3.h"

#include <SDL.h>

#include <atomic>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <queue>
#include <functional>
#include <condition_variable>

static std::function<bool()> g_doInit;
static std::function<bool()> g_doReload;
static std::function<void(int, int)> g_setWindowSize;
static std::function<int()> g_getOutputRecordId;
static std::function<bool()> g_mainUpdate;

void mainUpdate(void *) {
    g_mainUpdate();
}

// JS interface

extern "C" {
    EMSCRIPTEN_KEEPALIVE
        int do_init() {
            return g_doInit();
        }

    EMSCRIPTEN_KEEPALIVE
        int do_reload() {
            return g_doReload();
        }

    EMSCRIPTEN_KEEPALIVE
        void set_window_size(int sizeX, int sizeY) {
            g_setWindowSize(sizeX, sizeY);
        }

    EMSCRIPTEN_KEEPALIVE
        int get_output_id() {
            return g_getOutputRecordId();
        }
}

#ifdef __EMSCRIPTEN__
    int g_windowSizeX = 1400;
    int g_windowSizeY = 900;
#else
    int g_windowSizeX = 1920;
    int g_windowSizeY = 1200;
#endif

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

    int32_t keyPressWidth_samples   = kKeyWidth_samples;
    int32_t offsetFromPeak_samples  = kKeyWidth_samples - kKeyOffset_samples;
    int32_t alignWindow_samples     = kKeyAlign_samples;

    float wEnglish = 30.0f;
    int nHypothesesToKeep = 500;

    std::vector<int>   valuesClusters = { 40, 60, 80, 100, };
    std::vector<float> valuesFSpread  = { 0.5f, 0.75f, 1.0f, 1.5f, 2.0f, };

    int32_t nProcessors() const {
        return valuesClusters.size()*valuesFSpread.size();
    }

    int32_t valueForProcessorClusters(int idx) const {
      return valuesClusters[idx/(valuesFSpread.size())];
    }

    float valueForProcessorFSpread(int idx) const {
        return valuesFSpread[idx%valuesFSpread.size()];
    }

    Cipher::TParameters cipher;
};

struct stStateUI {
    bool doUpdate = false;

    struct Flags {
        bool recalculateSimilarityMap = false;
        bool resetOptimization = false;
        bool changeProcessing = false;
        bool applyParameters = false;
        bool applyHints = false;

        void clear() { memset(this, 0, sizeof(Flags)); }
    } flags;

    bool autoHint = false;
    bool processing = true;
    bool recording = false;
    bool audioCapture = false;
    bool calculatingSimilarityMap = false;

    // key presses window
    bool scrolling = false;
    bool recalculateKeyPresses = false;

    int nview = -1;
    int offset = -1;
    int viewMin = 512;
    int viewMax = 512;
    int nviewPrev = -1;
    int lastSize = -1;
    int lastKeyPresses = 0;

    float amin = std::numeric_limits<TSample>::min();
    float amax = std::numeric_limits<TSample>::max();

    float dragOffset = 0.0f;
    float scrollSize = 18.0f;

    TWaveform waveformLowRes;
    TWaveform waveformThreshold;
    TWaveform waveformMax;

    // window sizes
    float windowHeightTitleBar = 0;
    float windowHeightKeyPesses = 0;
    float windowHeightResults = 0;
    float windowHeightSimilarity = 0;

    bool openAboutWindow = false;
    bool openParametersWindow = false;
    bool loadRecord = false;
    bool loadKeyPresses = false;
    bool findBestCutoffFreq = false;
    bool rescaleWaveform = true;

    int outputRecordId = 0;

    double maxSample = 0.0;

    std::string fnameRecord = "default.kbd";
    std::string fnameKeyPressess = "default.kbd.keys";

    TParameters params;
    TWaveformF waveformOriginal;
    TWaveform waveformInput;
    TKeyPressCollection keyPresses;
    TSimilarityMap similarityMap;

    Cipher::THint suggestions;

    struct ProcessorResults : std::vector<Cipher::TResult> {
        int id = 0;
    };

    std::map<int, ProcessorResults> results;
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

struct stStateCapture {
    struct Flags {
        bool updateRecord = false;

        void clear() { memset(this, 0, sizeof(Flags)); }
    } flags;

    AudioLogger::Record recordNew;
};

template<typename T>
struct TripleBuffer : public T {
    bool update(bool force = false) {
        std::lock_guard<std::mutex> lock(mutex);
        hasChanged = true || force;
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

template<> bool TripleBuffer<stStateUI>::update(bool force) {
    std::lock_guard<std::mutex> lock(mutex);
    if (hasChanged && force == false) return false;

    hasChanged = true;

    buffer.flags = this->flags;

    if (this->flags.recalculateSimilarityMap || this->flags.resetOptimization) {
        buffer.params = this->params;
        buffer.keyPresses = this->keyPresses;
    }

    if (this->flags.changeProcessing) {
        buffer.processing = this->processing;
    }

    if (this->flags.applyParameters) {
        buffer.params = this->params;
    }

    if (this->flags.applyHints) {
        buffer.keyPresses = this->keyPresses;
    }

    this->flags.clear();

    return true;
}

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

template<> const stStateCapture & TripleBuffer<stStateCapture>::get() {
    std::lock_guard<std::mutex> lock(mutex);
    if (hasChanged) {
        out = std::move(buffer);
        hasChanged = false;
    }

    return out;
}

template<> bool TripleBuffer<stStateCapture>::update(bool force) {
    std::lock_guard<std::mutex> lock(mutex);
    if (hasChanged && force == false) return false;

    hasChanged = true;

    buffer.flags = this->flags;

    if (this->flags.updateRecord) {
        buffer.recordNew.insert(
                buffer.recordNew.end(),
                make_move_iterator(this->recordNew.begin()),
                make_move_iterator(this->recordNew.end()));
        this->recordNew.clear();
    }

    this->flags.clear();

    return true;
}

std::atomic<bool> g_doRecord;
AudioLogger audioLogger;
TripleBuffer<stStateUI> stateUI;
TripleBuffer<stStateCore> stateCore;
TripleBuffer<stStateCapture> stateCapture;

float plotWaveform(void * data, int i) {
    TWaveformView * waveform = (TWaveformView *)data;
    return waveform->samples[i];
}

float plotWaveformInverse(void * data, int i) {
    TWaveformView * waveform = (TWaveformView *)data;
    return -waveform->samples[i];
}

bool renderKeyPresses(stStateUI & stateUI, const TWaveform & waveform, TKeyPressCollection & keyPresses) {
    auto & params = stateUI.params;

    bool & scrolling = stateUI.scrolling;
    bool & recalculateKeyPresses = stateUI.recalculateKeyPresses;

    int & nview = stateUI.nview;
    int & viewMin = stateUI.viewMin;
    int & viewMax = stateUI.viewMax;
    int & offset = stateUI.offset;
    int & lastSize = stateUI.lastSize;
    int & lastKeyPresses = stateUI.lastKeyPresses;

    float & amin = stateUI.amin;
    float & amax = stateUI.amax;

    float & dragOffset = stateUI.dragOffset;
    float & scrollSize = stateUI.scrollSize;

    auto & nviewPrev = stateUI.nviewPrev;

    TWaveform & waveformLowRes = stateUI.waveformLowRes;
    TWaveform & waveformThreshold = stateUI.waveformThreshold;
    TWaveform & waveformMax = stateUI.waveformMax;

    if (nview < 0) nview = waveform.size();
    if (offset < 0) offset = (waveform.size() - nview)/2;

    if (stateUI.recording) {
        if (nview < 1024*kSamplesPerFrame && scrolling == false) {
            nview = std::min((int) (1024*kSamplesPerFrame), (int) waveform.size());
        }
    }

    if (lastSize != (int) waveform.size()) {
        viewMax = waveform.size();
        lastSize = waveform.size();
        nviewPrev = nview + 1;

        if (scrolling == false) {
            offset = waveform.size() - nview;
        }
    }

#if 0
    if ((int) keyPresses.size() >= lastKeyPresses + 10) {
        lastKeyPresses = keyPresses.size();

        stateUI.similarityMap.clear();
        stateUI.flags.recalculateSimilarityMap = true;
        stateUI.doUpdate = true;
    }
#endif

    ImGui::SetNextWindowPos(ImVec2(0, stateUI.windowHeightTitleBar), ImGuiCond_Once);
    if (stateUI.windowHeightKeyPesses == 0.0f) {
        ImGui::SetNextWindowSize(ImVec2(g_windowSizeX, 250.0f), ImGuiCond_Always);
    } else {
        ImGui::SetNextWindowSize(ImVec2(g_windowSizeX, stateUI.windowHeightKeyPesses), ImGuiCond_Always);
    }
    if (ImGui::Begin("Key Presses", nullptr, ImGuiWindowFlags_NoMove)) {
        bool ignoreDelete = false;

        auto wsize = ImVec2(ImGui::GetContentRegionAvailWidth(), ImGui::GetContentRegionAvail().y - 3*ImGui::GetTextLineHeightWithSpacing());

        if (nview != nviewPrev) {
            generateLowResWaveform(waveform, waveformLowRes, std::max(1.0f, nview/wsize.x));
            nviewPrev = nview;
        }

        auto wview = getView(waveformLowRes, offset, nview);

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

        if (waveform.size() == waveformThreshold.size() && waveform.size() == waveformMax.size()) {
            auto tview = getView(waveformThreshold, offset, nview);
            auto mview = getView(waveformMax, offset, nview);

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
        }
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

            if ((ImGui::IsMouseReleased(0) && ImGui::GetIO().KeyCtrl) || ImGui::IsMouseDoubleClicked(0)) {
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
                    if (((ImGui::IsMouseReleased(0) && ImGui::GetIO().KeyCtrl) || ImGui::IsMouseDoubleClicked(0)) && ignoreDelete == false) {
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

        static bool playHalfSpeed = false;
        static int historySizeReset = kFindKeysHistorySizeReset;
        static float thresholdBackground = kFindKeysThreshold;
        ImGui::PushItemWidth(100.0);

        if (stateUI.recording == false) {
            if (stateUI.audioCapture) {
                if (ImGui::Button("   Record")) {
                    if (float(stateUI.waveformInput.size()/kSampleRate) < kMaxRecordSize_s) {
                        audioLogger.resume();
                        stateUI.recording = true;
                        g_doRecord = true;
                    }
                }
                {
                    auto savePos = ImGui::GetCursorScreenPos();
                    auto p = savePos;
                    p.x += 0.90f*ImGui::CalcTextSize("  ").x;
                    p.y -= 0.56f*ImGui::GetFrameHeightWithSpacing();
                    ImGui::GetWindowDrawList()->AddCircleFilled(p, 6.0, ImGui::ColorConvertFloat4ToU32({ 1.0f, 0.0f, 0.0f, 1.0f }));
                    ImGui::SetCursorScreenPos(savePos);
                }

                ImGui::SameLine();
            }

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

            if (g_playbackData.idx >= g_playbackData.waveform.n - kSamplesPerFrame) {
                g_playbackData.playing = false;
                SDL_ClearQueuedAudio(g_deviceIdOut);
#ifndef __EMSCRIPTEN__
                SDL_PauseAudioDevice(g_deviceIdOut, 1);
#endif
            }

            ImGui::SameLine();
            ImGui::SliderFloat("Threshold background", &thresholdBackground, 0.1f, 50.0f) && (recalculateKeyPresses = true);
            ImGui::SameLine();
            ImGui::SliderInt("History Size", &historySizeReset, 512, 1024*16) && (recalculateKeyPresses = true);
            ImGui::SameLine();
            if (ImGui::Button("Recalculate")) {
                recalculateKeyPresses = true;
            }

            ImGui::SameLine();
            ImGui::DragInt("Key width", &params.keyPressWidth_samples, 8, 0, kSampleRate/10);
            ImGui::SameLine();
            ImGui::DragInt("Peak offset", &params.offsetFromPeak_samples, 8, -kSampleRate/10, kSampleRate/10);
            ImGui::SameLine();
            ImGui::DragInt("Align window", &params.alignWindow_samples, 8, 0, kSampleRate/10);

            ImGui::PopItemWidth();
        } else {
            if (g_doRecord) {
                g_doRecord = false;
                audioLogger.record(0.2f, 0);
            }

            if (ImGui::Button("   Stop Recording")) {
                audioLogger.pause();
                g_doRecord = false;

                stateUI.recalculateKeyPresses = true;
                stateUI.lastSize = stateUI.lastSize - 1;
                stateUI.recording = false;
                stateUI.findBestCutoffFreq = true;
                stateUI.rescaleWaveform = true;
            }
            {
                auto savePos = ImGui::GetCursorScreenPos();
                auto p0 = savePos;
                p0.x += 0.90f*ImGui::CalcTextSize("  ").x;
                p0.y -= 0.56f*ImGui::GetFrameHeightWithSpacing();

                auto p1 = p0;
                p0.x -= 4.0f;
                p0.y -= 4.0f;
                p1.x += 4.0f;
                p1.y += 4.0f;

                ImGui::GetWindowDrawList()->AddRectFilled(p0, p1, ImGui::ColorConvertFloat4ToU32({ 1.0f, 1.0f, 1.0f, 1.0f }));
                ImGui::SetCursorScreenPos(savePos);
            }

            ImGui::SameLine();
            ImGui::Text("Record length: %4.2f s (max : %4.2f s)", (float)(stateUI.waveformInput.size())/kSampleRate, kMaxRecordSize_s);
        }

        if (recalculateKeyPresses) {
            findKeyPresses(getView(waveform, 0), keyPresses, waveformThreshold, waveformMax, thresholdBackground, kFindKeysHistorySize, historySizeReset, kFindKeysRemoveLowPower);
            recalculateKeyPresses = false;
        }

        stateUI.windowHeightKeyPesses = ImGui::GetWindowHeight();
    } else {
        stateUI.windowHeightKeyPesses = ImGui::GetTextLineHeightWithSpacing();
    }
    ImGui::End();

    return false;
}

bool renderResults(stStateUI & stateUI) {
    float curHeight = std::min(g_windowSizeY - stateUI.windowHeightTitleBar - stateUI.windowHeightKeyPesses, (12 + std::min(8lu, stateUI.results.size())*kTopResultsPerProcessor)*ImGui::GetTextLineHeightWithSpacing());
    ImGui::SetNextWindowPos(ImVec2(0, stateUI.windowHeightTitleBar + stateUI.windowHeightKeyPesses), ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(1.0f*g_windowSizeX, curHeight), ImGuiCond_Always);
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
                stateUI.results.clear();
                if (stateUI.autoHint) {
                    for (auto & keyPress : stateUI.keyPresses) {
                        keyPress.bind = -1;
                    }
                }
                stateUI.flags.resetOptimization = true;
                stateUI.doUpdate = true;
            }

            //ImGui::SameLine();
            //if (ImGui::Button("Apply Suggestions")) {
            //    int n = stateUI.suggestions.size();
            //    for (int i = 0; i < n; ++i) {
            //        if (stateUI.suggestions[i] < 0) {
            //            stateUI.keyPresses[i].bind = -1;
            //            continue;
            //        }
            //        if (stateUI.suggestions[i] > 0 && stateUI.suggestions[i] <= 26) {
            //            stateUI.keyPresses[i].bind = stateUI.suggestions[i] - 1;
            //        } else {
            //            stateUI.keyPresses[i].bind = 26;
            //        }
            //    }

            //    stateUI.flags.applyHints = true;
            //    stateUI.doUpdate = true;
            //}

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

        //if (ImGui::SliderInt("Clusters", &stateUI.params.valuesClusters[0], 30, 150.0)) {
        //    stateUI.flags.applyParameters = true;
        //    stateUI.doUpdate = true;
        //}

        //if (ImGui::SliderFloat("Spread", &stateUI.params.fSpread, 0.1f, 3.0f)) {
        //    stateUI.flags.applyParameters = true;
        //    stateUI.doUpdate = true;
        //}

        ImGui::SameLine();
        if (ImGui::SliderInt("Keep", &stateUI.params.nHypothesesToKeep, 20, 2000)) {
            stateUI.flags.applyParameters = true;
            stateUI.doUpdate = true;
        }

        //ImGui::SameLine();
        //if (ImGui::Checkbox("Auto", &stateUI.autoHint)) {
        //}
        //if (ImGui::IsItemHovered()) {
        //    ImGui::BeginTooltip();
        //    ImGui::Text("If checked, hints will be updated automatically based on the suggestions");
        //    ImGui::EndTooltip();
        //}

        if (stateUI.results.size() > 0) {
            int nPerAutoHint = 10;
            auto id0 = stateUI.results.begin()->second.id%stateUI.params.cipher.nIters;
            auto id1 = (stateUI.results.begin()->second.id/stateUI.params.cipher.nIters)%nPerAutoHint;
            if (stateUI.processing) stateUI.results.begin()->second.id++;
            ImGui::SameLine();
            auto progress = float(id1)/nPerAutoHint;
            ImGui::ProgressBar(progress, { 100.0, ImGui::GetTextLineHeightWithSpacing() });
            if (stateUI.processing && id0 == 0 && id1 == 0) {
                // Important!
                // This decays the goodnes of found solutions in order to observe new solutions
                for (auto & result : stateUI.results) {
                    for (auto & item : result.second) {
                        item.p *= 1.01;
                    }
                }
                if (stateUI.autoHint && stateUI.results.begin()->second.id > 10000) {
                    int n = stateUI.suggestions.size();
                    for (int i = 0; i < n; ++i) {
                        if (frand() > 0.2) {
                            stateUI.keyPresses[i].bind = -1;
                            continue;
                        }
                        if (stateUI.suggestions[i] < 0) {
                            stateUI.keyPresses[i].bind = -1;
                            continue;
                        }

                        //if (stateUI.suggestions[i] == 0 || stateUI.suggestions[i] == 5) {
                        //    stateUI.keyPresses[i].bind = rand()%2 == 0 ? 26 : 4;
                        //    continue;
                        //}

                        //if (stateUI.suggestions[i] == 0 || stateUI.suggestions[i] == 5) {
                        //    stateUI.keyPresses[i].bind = -1;
                        //    continue;
                        //}

                        if (frand() > 0.0) {
                            if (stateUI.suggestions[i] > 0 && stateUI.suggestions[i] <= 26) {
                                stateUI.keyPresses[i].bind = stateUI.suggestions[i] - 1;
                            } else {
                                stateUI.keyPresses[i].bind = 26;
                            }
                        } else {
                            char c0 = (stateUI.suggestions[i] > 0 && stateUI.suggestions[i] <= 26) ?
                                'a' + stateUI.suggestions[i] - 1 : '_';

                            int idx = rand()%kNearbyKeys.at(c0).size();
                            char c1 = kNearbyKeys.at(c0)[idx];

                            if (c1 == '_') {
                                stateUI.keyPresses[i].bind = 26;
                            } else {
                                stateUI.keyPresses[i].bind = c1 - 'a';
                            }
                        }
                    }
                    stateUI.flags.applyHints = true;
                    stateUI.doUpdate = true;
                }
            }
        }

        ImGui::PopItemWidth();

        ImGui::BeginChild("##currentResults", { 0, 0 }, 1, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_AlwaysHorizontalScrollbar);

        ImGui::Separator();
        ImGui::Text("          %8s", "Suggestions");
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
                            c = 'a' + stateUI.suggestions[i] - 1;
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

        ImGui::Text("             %8s", "Hints");
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
                            c = 'a' + stateUI.keyPresses[i].bind;
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

                        if (ImGui::IsMouseDown(2)) {
                            stateUI.keyPresses[i].bind = 26;
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
                        ImGui::Text(" - Middle click to put a space as hint");
                        ImGui::EndTooltip();
                    }
                    if (i < n - 1) {
                        curPos.x += charSize.x;
                        ImGui::SetCursorScreenPos(curPos);
                    }
                }
            }
        }

        ImGui::Text("%2s. %8s %8s", "#", "Prob", "Clust");
        ImGui::Separator();

        for (const auto & items : stateUI.results) {
            ImGui::Separator();
            const auto & id = items.first;
            for (const auto & item : items.second) {
                const auto & result = item;

                ImGui::PushID(id);
                ImGui::Text("%2d. %8.3f %8.3f", id, result.p, result.pClusters);
                if (ImGui::IsItemHovered()) {
                    ImGui::BeginTooltip();
                    ImGui::Text("Processor   %d", id);
                    ImGui::Text("Cluster-to-letter map:");
                    {
                        int i = 0;
                        for (auto& pair : result.clMap) {
                            char c = '_';
                            if (pair.second > 0 && pair.second <= 26) {
                                c = 'a' + pair.second - 1;
                            }
                            ImGui::Text("%3d: %c  ", pair.first, c);
                            if (++i % 10 != 0) ImGui::SameLine();
                        }
                    }
                    ImGui::EndTooltip();
                }
                if (result.clMap.empty() == false) {
                    ImGui::SameLine();
                    int n = std::min(result.clusters.size(), stateUI.keyPresses.size());
                    auto charSize = ImGui::CalcTextSize("a");
                    auto curPos = ImGui::GetCursorScreenPos();
                    if (n > 0) {
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
                                    c = 'a' + let - 1;
                                } else if (let == 27) {
                                    c = '_';
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
                                    if (ImGui::IsMouseDown(0) && ImGui::GetIO().KeyCtrl) {
                                        for (int j = 0; j < n; ++j) {
                                            auto let2 = result.clMap.at(result.clusters[j]);
                                            if (let2 > 0 && let2 <= 26) {
                                                if (stateUI.keyPresses[j].bind == 26) {
                                                    stateUI.keyPresses[j].bind = -1;
                                                }
                                            } else {
                                                stateUI.keyPresses[j].bind = 26;
                                            }
                                        }
                                        stateUI.flags.applyHints = true;
                                        stateUI.doUpdate = true;
                                    } else if (ImGui::IsMouseDown(0)) {
                                        if (let > 0 && let <= 26) {
                                            stateUI.keyPresses[i].bind = let - 1;
                                        } else {
                                            stateUI.keyPresses[i].bind = 26;
                                        }
                                        stateUI.flags.applyHints = true;
                                        stateUI.doUpdate = true;
                                    } else if (ImGui::IsMouseDown(1)) {
                                        stateUI.keyPresses[i].bind = -1;
                                    }
                                    ImGui::BeginTooltip();
                                    ImGui::Text("Key presss: %d, letter - %d", i, let);
                                    ImGui::Text(" - Left click to set as hint");
                                    ImGui::Text(" - Ctrl + left click to set only the spaces as hint");
                                    ImGui::EndTooltip();
                                }
                            }

                            if (i < n - 1) {
                                curPos.x += charSize.x;
                                ImGui::SetCursorScreenPos(curPos);
                            }
                        }
                    }
                }
                ImGui::PopID();
            }
        }

        ImGui::EndChild();

        stateUI.windowHeightResults = ImGui::GetWindowHeight();
    } else {
        stateUI.windowHeightResults = ImGui::GetTextLineHeightWithSpacing();
    }

    ImGui::End();
    return true;
}

bool renderSimilarity(const TKeyPressCollection & keyPresses, const TSimilarityMap & similarityMap) {
    int offsetY = stateUI.windowHeightTitleBar + stateUI.windowHeightKeyPesses + stateUI.windowHeightResults;
    ImGui::SetNextWindowPos(ImVec2(0, offsetY), ImGuiCond_Always);
    ImGui::SetNextWindowSize(ImVec2(1.0f*g_windowSizeX, g_windowSizeY - offsetY), ImGuiCond_Always);
    if (ImGui::Begin("Similarity", nullptr, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_NoMove)) {
        auto wsize = ImGui::GetContentRegionAvail();

        static float bsize = 10.0f;
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
                            //ImGui::Text("[%3d, %3d]\n", keyPresses[i].cid, keyPresses[j].cid);
                            ImGui::Text("[%3d, %3d] = %5.4g\n", i, j, similarityMap[i][j].cc);
                            //for (int k = 0; k < n; ++k) {
                            //    if (similarityMap[i][k].cc > 0.5) ImGui::Text("Offset [%3d, %3d] = %d\n", i, k, (int) similarityMap[i][k].offset);
                            //}
                            //ImGui::Separator();
                            //for (int k = 0; k < n; ++k) {
                            //    if (similarityMap[k][i].cc > 0.5) ImGui::Text("Offset [%3d, %3d] = %d\n", k, i, (int) similarityMap[k][i].offset);
                            //}
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
        } else {
            ImGui::TextColored({1.0f, 0.3f, 0.4f, 1.0f}, "Waiting for similarity map to be available!");
        }

        ImGui::EndChild();

        stateUI.windowHeightSimilarity = ImGui::GetWindowHeight();
    } else {
        stateUI.windowHeightSimilarity = ImGui::GetTextLineHeightWithSpacing();
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

    printf("Build: %s, (%s)\n", kGIT_DATE, kGIT_SHA1);
    printf("Usage: %s record.kbd n-gram-dir [-pN] [-cN] [-CN] [-FN] [-fN]\n", argv[0]);
    printf("    -pN - select playback device N\n");
    printf("    -cN - select capture device N\n");
    printf("    -CN - select number N of capture channels to use\n");
    printf("    -FN - select filter type, (0 - none, 1 - first order high-pass, 2 - second order high-pass)\n");
    printf("    -fN - cutoff frequency in Hz\n");

    if (argc < 3) {
        return -1;
    }

    const auto argm = parseCmdArguments(argc, argv);
    const int playbackId    = argm.count("p") == 0 ? 0 : std::stoi(argm.at("p"));
    const int captureId     = argm.count("c") == 0 ? 0 : std::stoi(argm.at("c"));
    const int nChannels     = argm.count("C") == 0 ? 0 : std::stoi(argm.at("C"));
    const int filterId      = argm.count("F") == 0 ? EAudioFilter::FirstOrderHighPass : std::stoi(argm.at("F"));
    const int freqCutoff_Hz = argm.count("f") == 0 ? 0 : std::stoi(argm.at("f"));

    stateUI.params.playbackId = playbackId;
    stateUI.fnameRecord = argv[1];
    stateUI.fnameKeyPressess = stateUI.fnameRecord + ".keys";

    stateUI.waveformInput.reserve(kSamplesPerFrame*kMaxRecordSize_frames);
    stateUI.waveformOriginal.reserve(kSamplesPerFrame*kMaxRecordSize_frames);

    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) != 0) {
        printf("Error: %s\n", SDL_GetError());
        return -2;
    }

    AudioLogger::Callback cbAudio = [&](const AudioLogger::Record & frames) {
        stateCapture.recordNew.insert(
                stateCapture.recordNew.end(),
                frames.begin(),
                frames.end());

        stateCapture.flags.updateRecord = true;
        stateCapture.update();

        g_doRecord = true;

        return;
    };

    g_doInit = [&]() {
        if (prepareAudioOut(stateUI.params) == false) {
            printf("Error: failed to initialize audio playback\n");
            return false;
        }

        AudioLogger::Parameters parameters;
        parameters.callback = cbAudio;
        parameters.captureId = captureId;
        parameters.nChannels = nChannels;
        parameters.sampleRate = kSampleRate;
        parameters.filter = EAudioFilter::None;
        parameters.freqCutoff_Hz = kFreqCutoff_Hz;

        if (audioLogger.install(std::move(parameters)) == false) {
            fprintf(stderr, "Failed to initialize audio capture\n");
        } else {
            audioLogger.pause();
            stateUI.audioCapture = true;

            printf("[+] Audio capture initilized succesfully\n");
        }

        return true;
    };

    g_doReload = [&]() {
        stateUI.loadRecord = true;

        return true;
    };

    g_getOutputRecordId = [&]() {
        return stateUI.outputRecordId;
    };

    printf("[+] Loading recording from '%s'\n", argv[1]);
    if (readFromFile<TSampleF>(argv[1], stateUI.waveformOriginal) == false) {
        printf("Specified file '%s' does not exist\n", argv[1]);
        //return -4;
    } else {
        auto freqCutoffCur_Hz = freqCutoff_Hz;
        if (freqCutoffCur_Hz == 0) {
            const auto tStart = std::chrono::high_resolution_clock::now();

            freqCutoffCur_Hz = Cipher::findBestCutoffFreq(stateUI.waveformOriginal, (EAudioFilter) filterId, kSampleRate, 100.0f, 1000.0f, 100.0f);

            const auto tEnd = std::chrono::high_resolution_clock::now();
            printf("[+] Found best freqCutoff = %d Hz, took %4.3f seconds\n", freqCutoffCur_Hz, toSeconds(tStart, tEnd));
        }

        printf("[+] Filtering waveform with filter type = %d and cutoff frequency = %d Hz\n", filterId, freqCutoffCur_Hz);
        ::filter(stateUI.waveformOriginal, (EAudioFilter) filterId, freqCutoffCur_Hz, kSampleRate);

        printf("[+] Converting waveform to i16 format ...\n");
        if (convert(stateUI.waveformOriginal, stateUI.waveformInput) == false) {
            printf("Conversion failed\n");
            return -4;
        }

        stateUI.recalculateKeyPresses = true;
        stateUI.maxSample = calcAbsMax(stateUI.waveformOriginal);
    }

    Cipher::TFreqMap freqMap6;

    if (Cipher::loadFreqMapBinary((std::string(argv[2]) + "/./ggwords-6-gram.dat.binary").c_str(), freqMap6) == false) {
        return -5;
    }

    stateCore.freqMap[0] = &freqMap6;
    stateCore.freqMap[1] = &freqMap6;
    stateCore.freqMap[2] = &freqMap6;

    Gui::Objects guiObjects;
    if (Gui::init("Keytap3", g_windowSizeX, g_windowSizeY, guiObjects) == false) {
        return -6;
    }

    g_setWindowSize = [&](int sizeX, int sizeY) {
        SDL_SetWindowSize(guiObjects.window, sizeX, sizeY);
    };

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
                    if (stateCoreNew.processors[i].getResult().id != stateUI.results[i].id) {
                        stateUI.results[i].id = stateCoreNew.processors[i].getResult().id;
                        if (stateUI.results[i].size() < kTopResultsPerProcessor) {
                            stateUI.results[i].push_back(stateCoreNew.processors[i].getResult());
                        } else if (stateCoreNew.processors[i].getResult().p > stateUI.results[i].back().p) {
                            stateUI.results[i].back() = stateCoreNew.processors[i].getResult();
                        }

                        int k = stateUI.results[i].size() - 1;
                        while (k > 0) {
                            if (stateUI.results[i][k-1].p < stateUI.results[i][k].p) {
                                std::swap(stateUI.results[i][k-1], stateUI.results[i][k]);
                            }
                            --k;
                        }
                    }
                }
            }

            recalcSuggestions = false;
            if (recalcSuggestions) {
                int n = stateUI.keyPresses.size();
                stateUI.suggestions.resize(n);

                for (int i = 0; i < n; ++i) {
                    int c = 0;
                    std::map<int, int> cnt;
                    for (const auto & result : stateUI.results) {
                        const auto & item = result.second.front();
                        if (n != (int) item.clusters.size()) continue;
                        if (item.clMap.find(item.clusters[i]) == item.clMap.end()) continue;
                        if (i <= (int) item.clusters.size()) {
                            ++cnt[item.clMap.at(item.clusters[i])];
                            ++c;
                        }
                    }
                    bool hasSuggestion = false;
                    for (auto & s : cnt) {
                        if (s.second > 0.5*c) {
                            stateUI.suggestions[i] = s.first;
                            hasSuggestion = true;
                        }
                    }
                    if (hasSuggestion == false) {
                        stateUI.suggestions[i] = -1;
                    }
                }
            } else {
                int n = stateUI.keyPresses.size();
                stateUI.suggestions.resize(n);
            }
        }

        if (stateCapture.changed()) {
            auto stateCaptureNew = stateCapture.get();

            if (stateCaptureNew.flags.updateRecord) {
                for (auto & frame : stateCaptureNew.recordNew) {
                    stateUI.waveformOriginal.insert(stateUI.waveformOriginal.end(), frame.begin(), frame.end());
                }

                int nold = stateUI.waveformInput.size();
                float scale = float(std::numeric_limits<TSample>::max())/(stateUI.maxSample == 0.0 ? 1.0f : stateUI.maxSample);
                stateUI.waveformInput.resize(stateUI.waveformOriginal.size());
                for (int i = nold; i < (int) stateUI.waveformInput.size(); ++i) {
                    stateUI.waveformInput[i] = scale*stateUI.waveformOriginal[i];
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
                case SDL_DROPFILE:
                    {
                        char * path = event.drop.file;
                        stateUI.fnameRecord = path;
                        stateUI.loadRecord = true;
                        SDL_free(path);
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

        stateUI.windowHeightTitleBar = ImGui::GetTextLineHeightWithSpacing();

        if (ImGui::BeginMainMenuBar())
        {
            if (ImGui::BeginMenu("File"))
            {
                ImGui::Separator();
                ImGui::TextDisabled("Record: %s", stateUI.fnameRecord.c_str());
                ImGui::Separator();
                if (ImGui::MenuItem("Save Record")) {
                    saveToFile(stateUI.fnameRecord, stateUI.waveformOriginal);
                    stateUI.outputRecordId++;
                }
                if (ImGui::MenuItem("Load Record")) {
                    stateUI.loadRecord = true;
                }

                ImGui::Separator();
                ImGui::TextDisabled("Key Presses: %s", stateUI.fnameKeyPressess.c_str());
                ImGui::Separator();
                if (ImGui::MenuItem("Save Key Presses")) {
                    saveKeyPresses(stateUI.fnameKeyPressess.c_str(), stateUI.keyPresses);
                }
                if (ImGui::MenuItem("Load Key Presses")) {
                    stateUI.loadKeyPresses = true;
                }

                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Help")) {
                if (ImGui::MenuItem("Parameters")) {
                    stateUI.openParametersWindow = true;
                }
                if (ImGui::MenuItem("About")) {
                    stateUI.openAboutWindow = true;
                }
                ImGui::EndMenu();
            }
            ImGui::EndMainMenuBar();
        }

        if (stateUI.openAboutWindow) {
            ImGui::OpenPopup("About");
            stateUI.openAboutWindow = false;
        }
        ImGui::SetNextWindowSize({512, -1}, ImGuiCond_Once);
        if (ImGui::BeginPopupModal("About", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoScrollbar)) {
            ImGui::Text("%s", "");
            {
                const char * text = "Keytap3 : v0.1 [%s]    ";
                ImGui::Text("%s", "");
                ImGui::SameLine(0.5f*ImGui::GetContentRegionAvailWidth() - 0.45f*ImGui::CalcTextSize(text).x);
                ImGui::Text(text, kGIT_SHA1);
            }
            ImGui::Text("%s", "");
            {
                const char * text = "author : Georgi Gerganov";
                ImGui::Text("%s", "");
                ImGui::SameLine(0.5f*ImGui::GetContentRegionAvailWidth() - 0.45f*ImGui::CalcTextSize(text).x);
                ImGui::Text("%s", text);
            }
            {
                const char * text = "source code : https://github.com/ggerganov/kbd-audio";
                ImGui::Text("%s", "");
                ImGui::SameLine(0.5f*ImGui::GetContentRegionAvailWidth() - 0.45f*ImGui::CalcTextSize(text).x);
                ImGui::Text("%s", text);
            }
            ImGui::Text("%s", "");
            {
                const char * text = "Close";
                ImGui::Text("%s", "");
                ImGui::SameLine(0.5f*ImGui::GetContentRegionAvailWidth() - 0.45f*ImGui::CalcTextSize(text).x);
                if (ImGui::Button(text)) {
                    ImGui::CloseCurrentPopup();
                }
            }
            ImGui::Text("%s", "");
            ImGui::EndPopup();
        }

        if (stateUI.openParametersWindow) {
            ImGui::OpenPopup("Parameters");
            stateUI.openParametersWindow = false;
        }
        ImGui::SetNextWindowSize({320, -1}, ImGuiCond_Once);
        if (ImGui::BeginPopupModal("Parameters", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoScrollbar)) {
            ImGui::Text("%s", "");
            ImGui::Text("Sample rate            : %d\n", (int) kSampleRate);
            ImGui::Text("Samples per frame      : %d\n", (int) kSamplesPerFrame);
            ImGui::Text("Hardware concurrency   : %d\n", (int) std::thread::hardware_concurrency());
            ImGui::Text("Max record size        : %g s\n", kMaxRecordSize_s);
            ImGui::Text("Freq cutoff            : %g Hz\n", kFreqCutoff_Hz);
            ImGui::Text("Subbreak Processors    : %d\n", stateUI.params.nProcessors());
            ImGui::Text("%s", "");

            {
                const char * text = "Close";
                ImGui::Text("%s", "");
                ImGui::SameLine(0.5f*ImGui::GetContentRegionAvailWidth() - 0.45f*ImGui::CalcTextSize(text).x);
                if (ImGui::Button(text)) {
                    ImGui::CloseCurrentPopup();
                }
            }
            ImGui::Text("%s", "");
            ImGui::EndPopup();
        }


        if (stateUI.loadRecord) {
            printf("[+] Loading recording from '%s'\n", argv[1]);
            if (readFromFile<TSampleF>(stateUI.fnameRecord, stateUI.waveformOriginal) == false) {
                printf("Specified file '%s' does not exist\n", argv[1]);
            } else {
                auto freqCutoffCur_Hz = freqCutoff_Hz;
                if (freqCutoffCur_Hz == 0) {
                    const auto tStart = std::chrono::high_resolution_clock::now();

                    freqCutoffCur_Hz = Cipher::findBestCutoffFreq(stateUI.waveformOriginal, (EAudioFilter) filterId, kSampleRate, 100.0f, 1000.0f, 100.0f);

                    const auto tEnd = std::chrono::high_resolution_clock::now();
                    printf("[+] Found best freqCutoff = %d Hz, took %4.3f seconds\n", freqCutoffCur_Hz, toSeconds(tStart, tEnd));
                }

                printf("[+] Filtering waveform with filter type = %d and cutoff frequency = %d Hz\n", filterId, freqCutoffCur_Hz);
                ::filter(stateUI.waveformOriginal, (EAudioFilter) filterId, freqCutoffCur_Hz, kSampleRate);

                printf("[+] Converting waveform to i16 format ...\n");
                if (convert(stateUI.waveformOriginal, stateUI.waveformInput) == false) {
                    printf("Conversion failed\n");
                } else {
                    stateUI.nview = -1;
                    stateUI.recalculateKeyPresses = true;
                    stateUI.maxSample = calcAbsMax(stateUI.waveformOriginal);
                    stateUI.flags.recalculateSimilarityMap = true;
                    stateUI.doUpdate = true;
                }
            }

            stateUI.loadRecord = false;
        }

        if (stateUI.loadKeyPresses) {
            loadKeyPresses(stateUI.fnameKeyPressess.c_str(), getView(stateUI.waveformInput, 0), stateUI.keyPresses);
            stateUI.loadKeyPresses = false;
        }

        if (stateUI.findBestCutoffFreq) {
            auto freqCutoffCur_Hz = freqCutoff_Hz;
            if (freqCutoffCur_Hz == 0) {
                const auto tStart = std::chrono::high_resolution_clock::now();

                freqCutoffCur_Hz = Cipher::findBestCutoffFreq(stateUI.waveformOriginal, (EAudioFilter) filterId, kSampleRate, 100.0f, 1000.0f, 100.0f);

                const auto tEnd = std::chrono::high_resolution_clock::now();
                printf("[+] Found best freqCutoff = %d Hz, took %4.3f seconds\n", freqCutoffCur_Hz, toSeconds(tStart, tEnd));
            }

            printf("[+] Filtering waveform with filter type = %d and cutoff frequency = %d Hz\n", filterId, freqCutoffCur_Hz);
            ::filter(stateUI.waveformOriginal, (EAudioFilter) filterId, freqCutoffCur_Hz, kSampleRate);

            stateUI.findBestCutoffFreq = false;
        }

        if (stateUI.rescaleWaveform) {
            if (convert(stateUI.waveformOriginal, stateUI.waveformInput) == false) {
                fprintf(stderr, "error : recording failed\n");
            }
            stateUI.maxSample = calcAbsMax(stateUI.waveformOriginal);

            stateUI.rescaleWaveform = false;
        }

        if (float(stateUI.waveformInput.size())/kSampleRate >= kMaxRecordSize_s) {
            stateUI.waveformInput.resize((kMaxRecordSize_frames - 1)*kSamplesPerFrame);
            audioLogger.pause();
            g_doRecord = false;

            stateUI.recalculateKeyPresses = true;
            stateUI.lastSize = stateUI.lastSize - 1;
            stateUI.recording = false;
            stateUI.rescaleWaveform = true;
        }

        for (int i = 0; i < IM_ARRAYSIZE(ImGui::GetIO().KeysDown); i++) if (ImGui::IsKeyPressed(i)) {
            //printf("Key pressed: %d\n", i);
        }

        renderKeyPresses(stateUI, stateUI.waveformInput, stateUI.keyPresses);
        renderResults(stateUI);
        renderSimilarity(stateUI.keyPresses, stateUI.similarityMap);

        Gui::render(guiObjects);

        if (stateUI.doUpdate) {
            stateUI.update();
            stateUI.doUpdate = false;
        }

        return true;
    };

    std::thread workerCore([&]() {
        while (finishApp == false) {
            if (stateUI.changed()) {
                auto stateUINew = stateUI.get();

                if (stateUINew.flags.recalculateSimilarityMap || stateUINew.flags.resetOptimization) {
                    if (stateUINew.keyPresses.size() < 3) continue;

                    stateCore.params = stateUINew.params;
                    stateCore.keyPresses = stateUINew.keyPresses;

                    printf("[+] Recalculating similarity map ...\n");

                    if (stateUINew.flags.recalculateSimilarityMap) {
                        stateCore.flags.calculatingSimilarityMap = true;
                        stateCore.update(true);

                        calculateSimilartyMap(
                                stateUINew.params.keyPressWidth_samples,
                                stateUINew.params.alignWindow_samples,
                                stateUINew.params.offsetFromPeak_samples,
                                stateCore.keyPresses,
                                stateCore.similarityMap);

                        int nTries = 3;
                        while (adjustKeyPresses(stateCore.keyPresses, stateCore.similarityMap) && --nTries) {
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
                        stateCore.update(true);
                    }

                    int n = stateCore.params.nProcessors();
                    for (int i = 0; i < n; ++i) {
                        const auto nClusters = stateCore.params.valueForProcessorClusters(i);
                        const auto wEnglish = stateCore.params.wEnglish;
                        const auto fSpread = stateCore.params.valueForProcessorFSpread(i);
                        const auto nHypothesesToKeep = stateCore.params.nHypothesesToKeep;

                        Cipher::TParameters params;
                        params.maxClusters = nClusters;
                        params.wEnglishFreq = wEnglish;
                        params.fSpread = fSpread;
                        params.nHypothesesToKeep = nHypothesesToKeep;
                        stateCore.processors[i] = Cipher::Processor();
                        stateCore.processors[i].init(
                                params,
                                *stateCore.freqMap[i%3],
                                stateCore.similarityMap);

                        printf("[+] Processor %d initialized: cluster = %d, wEnglish = %g, fSpread = %g, nHypothesesToKeep = %d\n",
                                i, nClusters, wEnglish, fSpread, nHypothesesToKeep);
                    }
                }

                if (stateUINew.flags.changeProcessing) {
                    stateCore.processing = stateUINew.processing;
                }

                if (stateUINew.flags.applyParameters) {
                    stateCore.params = stateUINew.params;

                    int n = stateCore.params.nProcessors();
                    for (int i = 0; i < n; ++i) {
                        const auto nClusters = stateCore.params.valueForProcessorClusters(i);
                        const auto wEnglish = stateCore.params.valueForProcessorFSpread(i);
                        const auto fSpread = stateCore.params.valueForProcessorFSpread(i);
                        const auto nHypothesesToKeep = stateCore.params.nHypothesesToKeep;

                        Cipher::TParameters params;
                        params.maxClusters = nClusters;
                        params.wEnglishFreq = wEnglish;
                        params.fSpread = fSpread;
                        params.nHypothesesToKeep = nHypothesesToKeep;
                        stateCore.processors[i].init(
                                params,
                                *stateCore.freqMap[i%3],
                                stateCore.similarityMap);
                    }
                }

                if (stateUINew.flags.applyHints) {
                    stateCore.keyPresses = stateUINew.keyPresses;
                }
            }

            if (stateCore.processing && stateCore.keyPresses.size() >= 3) {
                {
                    int n = stateCore.keyPresses.size();
                    stateCore.params.cipher.hint.clear();
                    stateCore.params.cipher.hint.resize(n, -1);
                    for (int i = 0; i < n; ++i) {
                        if (stateCore.keyPresses[i].bind < 0) continue;
                        stateCore.params.cipher.hint[i] = stateCore.keyPresses[i].bind + 1;
                    }
                }

#ifdef __EMSCRIPTEN__
                static int iter = 0;
                {
                    int p = iter%stateCore.params.nProcessors();
                    stateCore.processors[p].setHint(stateCore.params.cipher.hint);
                    stateCore.processors[p].compute();

                    stateCore.flags.updateResult[p] = true;
                    stateCore.update();

                    ++iter;
                }
#else
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
#endif
            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
    });

#ifdef __EMSCRIPTEN__
    emscripten_set_main_loop_arg(mainUpdate, NULL, 0, true);
#else
    if (g_doInit() == false) {
        printf("Error: failed to initialize audio playback\n");
        return -2;
    }

    while (true) {
        if (g_mainUpdate() == false) break;
    }

    workerCore.join();

    printf("[+] Terminated\n");

    Gui::free(guiObjects);
#endif

    return 0;
}
