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
#include <thread>

#define MY_DEBUG

constexpr int64_t kKeyPressWidth_samples = 256;
constexpr int64_t kSampleRate = 24000;

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
    TCC     cc = 0.0;
    TOffset offset = 0;
};

struct stWaveformView {
    const TSample * samples = nullptr;
    int64_t         n = 0;
};

struct stKeyPressData {
    TWaveformView       waveform;
    TKeyPressPosition   pos = 0;
    TCC                 ccAvg = 0.0;
};

template <typename T>
float toSeconds(T t0, T t1) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()/1024.0f;
}

inline float frand() { return ((float)rand())/RAND_MAX; }

TWaveformView getView(const TWaveform & waveform, int64_t idx) { return { waveform.data() + idx, (int64_t) waveform.size() - idx }; }

TWaveformView getView(const TWaveform & waveform, int64_t idx, int64_t len) { return { waveform.data() + idx, len }; }

bool saveKeyPresses(const char * fname, const TKeyPressCollection & keyPresses) {
    std::ofstream fout(fname, std::ios::binary);
    int n = keyPresses.size();
    fout.write((char *)(&n), sizeof(n));
    for (int i = 0; i < n; ++i) {
        fout.write((char *)(&keyPresses[i].pos), sizeof(keyPresses[i].pos));
    }
    fout.close();

    return true;
}

bool loadKeyPresses(const char * fname, const TWaveformView & waveform, TKeyPressCollection & keyPresses) {
    keyPresses.clear();

    std::ifstream fin(fname, std::ios::binary);
    int n = 0;
    fin.read((char *)(&n), sizeof(n));
    keyPresses.resize(n);
    for (int i = 0; i < n; ++i) {
        keyPresses[i].waveform = waveform;
        fin.read((char *)(&keyPresses[i].pos), sizeof(keyPresses[i].pos));
    }
    fin.close();

    return true;
}

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

bool generateLowResWaveform(const TWaveformView & waveform, TWaveform & waveformLowRes, int nWindow) {
    waveformLowRes.resize(waveform.n);

    int k = nWindow;

    std::deque<int64_t> que(k);
    auto [samples, n] = waveform;

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

bool findKeyPresses(const TWaveformView & waveform, TKeyPressCollection & res, TWaveform & waveformThreshold, double thresholdBackground, int historySize) {
    res.clear();
    waveformThreshold.resize(waveform.n);

    int rbBegin = 0;
    double rbAverage = 0.0;
    std::vector<double> rbSamples(8*historySize, 0.0);

    int k = historySize;

    std::deque<int64_t> que(k);
    auto [samples, n] = waveform;

    TWaveform waveformAbs(n);
    for (int64_t i = 0; i < n; ++i) {
        waveformAbs[i] = std::abs(samples[i]);
    }

    for (int64_t i = 0; i < n; ++i) {
        {
            int64_t ii = i - k/2;
            if (ii >= 0) {
                rbAverage *= rbSamples.size();
                rbAverage -= rbSamples[rbBegin];
                double acur = waveformAbs[i];
                rbSamples[rbBegin] = acur;
                rbAverage += acur;
                rbAverage /= rbSamples.size();
                if (++rbBegin >= rbSamples.size()) {
                    rbBegin = 0;
                }
            }
        }

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
            if (itest >= 2*k && itest < n - 2*k && que.front() == itest) {
                double acur = waveformAbs[itest];
                if (acur > thresholdBackground*rbAverage){
                    TKeyPressData entry;
                    entry.waveform = waveform;
                    entry.pos = itest;
                    entry.ccAvg = 0.0;
                    res.emplace_back(std::move(entry));
                }
            }
            waveformThreshold[itest] = waveformAbs[que.front()];
        }
    }

    return true;
}

bool findKeyPresses(const TWaveform & waveform, TKeyPressCollection & res, TWaveform & waveformThreshold, double thresholdBackground = 10.0, int historySize = 4*1024) {
    return findKeyPresses(getView(waveform, 0), res, waveformThreshold, thresholdBackground, historySize);
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
    int64_t idx = 0;
    int64_t offset = 0;
    TWaveformView waveform;
};

SDL_AudioDeviceID g_deviceIdOut = 0;
PlaybackData g_playbackData;

bool renderKeyPresses(const char * fnameInput, const TWaveform & waveform, TKeyPressCollection & keyPresses) {
    if (ImGui::Begin("Key Presses")) {
        int viewMin = 512;
        int viewMax = waveform.size();

        bool ignoreDelete = false;

        static int nview = waveform.size()/16;
        static int offset = (waveform.size() - nview)/2;
        static float amin = -16000;
        static float amax = 16000;
        static float dragOffset = 0.0f;
        static float scrollSize = 18.0f;

        static auto nviewPrev = nview + 1;

        static TWaveform waveformLowRes = waveform;
        static TWaveform waveformThreshold = waveform;

        auto wview = getView(waveformLowRes, offset, nview);
        auto tview = getView(waveformThreshold, offset, nview);
        auto wsize = ImVec2(ImGui::GetContentRegionAvailWidth(), 250.0);

        auto mpos = ImGui::GetIO().MousePos;
        auto savePos = ImGui::GetCursorScreenPos();
        auto drawList = ImGui::GetWindowDrawList();
        ImGui::PushStyleColor(ImGuiCol_FrameBg, { 0.1f, 0.1f, 0.1f, 0.3f });
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
        ImGui::PushStyleColor(ImGuiCol_PlotLines, { 0.0f, 1.0f, 0.0f, 0.5f });
        ImGui::PlotLines("##WaveformThreshold", plotWaveform, &tview, nview, 0, "", amin, amax, wsize);
        ImGui::PopStyleColor(2);
        ImGui::SetCursorScreenPos(savePos);
        ImGui::PushStyleColor(ImGuiCol_FrameBg, { 0.0f, 0.0f, 0.0f, 0.0f });
        ImGui::PushStyleColor(ImGuiCol_PlotLines, { 0.0f, 1.0f, 0.0f, 0.5f });
        ImGui::PlotLines("##WaveformThreshold", plotWaveformInverse, &tview, nview, 0, "", amin, amax, wsize);
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
                int pos = offset + nview*(mpos.x - savePos.x)/wsize.x;
                for (i = 0; i < keyPresses.size(); ++i) {
                    if (std::abs(keyPresses[i].pos - pos) < kKeyPressWidth_samples) break;
                    if (keyPresses[i].pos > pos) {
                        ignoreDelete = true;
                        TKeyPressData entry;
                        entry.pos = pos;
                        entry.waveform = wview;
                        keyPresses.insert(keyPresses.begin() + i, entry);
                        break;
                    }
                }
                if (i == keyPresses.size() && ignoreDelete == false) {
                    ignoreDelete = true;
                    TKeyPressData entry;
                    entry.pos = pos;
                    entry.waveform = wview;
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
            if (keyPresses[i].pos + kKeyPressWidth_samples < offset) continue;
            if (keyPresses[i].pos - kKeyPressWidth_samples >= offset + nview) break;

            {
                float x0 = ((float)(keyPresses[i].pos - kKeyPressWidth_samples - offset))/nview;
                float x1 = ((float)(keyPresses[i].pos + kKeyPressWidth_samples - offset))/nview;
                float x2 = ((float)(g_playbackData.offset + g_playbackData.idx - offset))/nview;

                ImVec2 p0 = {savePos.x + x0*wsize.x, savePos.y};
                ImVec2 p1 = {savePos.x + x1*wsize.x, savePos.y + wsize.y};
                bool isHovered = (mpos.x > p0.x && mpos.x < p1.x && mpos.y > p0.y && mpos.y < p1.y);

                float col = isHovered ? 0.7f : 0.3f;
                drawList->AddRectFilled(p0, p1, ImGui::ColorConvertFloat4ToU32({1.0f, 0.0f, 0.0f, col}));
                drawList->AddRect({savePos.x + x2*wsize.x, savePos.y}, {savePos.x + x2*wsize.x + 1.0f, savePos.y + wsize.y}, ImGui::ColorConvertFloat4ToU32({1.0f, 1.0f, 0.0f, 0.3f}));

                if (isHovered) {
                    if (ImGui::IsMouseReleased(0) && ImGui::GetIO().KeyCtrl && ignoreDelete == false) {
                        keyPresses.erase(keyPresses.begin() + i);
                        --i;
                    }
                }

                if (nview < 64.0*wsize.x) {
                    ImGui::SetCursorScreenPos({savePos.x + 0.5f*((x0 + x1)*wsize.x - ImGui::CalcTextSize(std::to_string(i).c_str()).x), savePos.y + wsize.y - ImGui::GetTextLineHeightWithSpacing()});
                    ImGui::Text("%d", i);
                }
            }
        }
        for (int i = 0; i < (int) keyPresses.size(); ++i) {
            {
                float x0 = ((float)(keyPresses[i].pos))/viewMax;
                drawList->AddRect({savePos.x + x0*wsize.x, savePos.y + wsize.y}, {savePos.x + x0*wsize.x + 1.0f, savePos.y + wsize.y + scrollSize}, ImGui::ColorConvertFloat4ToU32({1.0f, 0.0f, 0.0f, 0.3f}));
            }
        }

        ImGui::SetCursorScreenPos(savePos2);

        if (ImGui::Button("Play")) {
            g_playbackData.idx = 0;
            g_playbackData.offset = offset;
            g_playbackData.waveform = getView(waveform, offset, std::min((int) (10*kSampleRate), nview));
            SDL_PauseAudioDevice(g_deviceIdOut, 0);
        }

        if (g_playbackData.idx > g_playbackData.waveform.n) {
            SDL_PauseAudioDevice(g_deviceIdOut, 1);
        }

        ImGui::SameLine();
        static bool recalculate = true;
        static int historySize = 6*1024;
        static float thresholdBackground = 10.0;
        ImGui::PushItemWidth(100.0);
        ImGui::SliderFloat("Threshold background", &thresholdBackground, 0.1f, 50.0f) && (recalculate = true);
        ImGui::SameLine();
        ImGui::SliderInt("History Size", &historySize, 512, 1024*16) && (recalculate = true);
        ImGui::SameLine();
        if (ImGui::Button("Recalculate") || recalculate) {
            findKeyPresses(waveform, keyPresses, waveformThreshold, thresholdBackground, historySize);
            recalculate = false;
        }
        ImGui::PopItemWidth();

        static std::string filename = std::string(fnameInput) + ".keys";
        ImGui::SameLine();
        ImGui::Text("%s", filename.c_str());

        ImGui::SameLine();
        if (ImGui::Button("Save")) {
            saveKeyPresses(filename.c_str(), keyPresses);
        }

        ImGui::SameLine();
        if (ImGui::Button("Load")) {
            loadKeyPresses(filename.c_str(), wview, keyPresses);
        }

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
    auto end = std::min(data->idx + PlaybackData::kSamples, data->waveform.n);
    auto idx = data->idx;
    for (; idx < end; ++idx) {
        int16_t a = data->waveform.samples[idx];
        memcpy(stream + (idx - data->idx)*sizeof(a), &a, sizeof(a));
        len -= sizeof(a);
    }
    while (len > 0) {
        int16_t a = 0;
        memcpy(stream + (idx - data->idx)*sizeof(a), &a, sizeof(a));
        len -= sizeof(a);
        ++idx;
    }
    data->idx = idx;
}

bool prepareAudioOut() {
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

    playbackSpec.freq = kSampleRate;
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

    printf("Opened playback device %d\n", g_deviceIdOut);
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

    if (SDL_Init(SDL_INIT_VIDEO|SDL_INIT_TIMER) != 0) {
        printf("Error: %s\n", SDL_GetError());
        return -1;
    }

    if (prepareAudioOut() == false) {
        printf("Error: failed to initialize audio playback\n");
        return -2;
    }

    int windowSizeX = 1920;
    int windowSizeY = 1024;

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
    SDL_Window* window = SDL_CreateWindow("Keytap", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, windowSizeX, windowSizeY, SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE|SDL_WINDOW_ALLOW_HIGHDPI);
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

    int64_t sampleRate = kSampleRate;

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
    //{
    //    auto tStart = std::chrono::high_resolution_clock::now();
    //    printf("[+] Searching for key presses\n");
    //    if (findKeyPresses(waveformInput, keyPresses) == false) {
    //        printf("Failed to detect keypresses\n");
    //        return -2;
    //    }
    //    auto tEnd = std::chrono::high_resolution_clock::now();
    //    printf("[+] Detected a total of %d potential key presses\n", (int) keyPresses.size());
    //    for (auto & k : keyPresses) {
    //        printf("    position - %d\n", (int) k.pos);
    //    }
    //    printf("[+] Search took %4.3f seconds\n", toSeconds(tStart, tEnd));

    //    dumpKeyPresses("key_presses.plot", keyPresses);
    //}

    TSimilarityMap similarityMap;
    //{
    //    auto tStart = std::chrono::high_resolution_clock::now();
    //    printf("[+] Calculating CC similarity map\n");
    //    if (calculateSimilartyMap(keyPresses, similarityMap) == false) {
    //        printf("Failed to calculate similariy map\n");
    //        return -3;
    //    }
    //    auto tEnd = std::chrono::high_resolution_clock::now();
    //    printf("[+] Calculation took %4.3f seconds\n", toSeconds(tStart, tEnd));
    //}

    int n = keyPresses.size();
    for (int i = 0; i < n; ++i) {
        auto avgcc = keyPresses[i].ccAvg;
        printf("    Average CC for keypress %4d - %6.3f\n", i, avgcc);
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

        auto tStart = std::chrono::high_resolution_clock::now();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplSDL2_NewFrame(window);
        ImGui::NewFrame();

        renderKeyPresses(argv[1], waveformInput, keyPresses);

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
