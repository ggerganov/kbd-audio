/*! \file common-gui.cpp
 *  \brief Enter description here.
 */

#include "common-gui.h"

#include <deque>

template<typename T>
bool generateLowResWaveform(const TWaveformViewT<T> & waveform, TWaveformT<T> & waveformLowRes, int nWindow) {
    waveformLowRes.resize(waveform.n);

    int k = nWindow;
    std::deque<int64_t> que(k);

    //auto [samples, n] = waveform;
    auto samples = waveform.samples;
    auto n       = waveform.n;

    TWaveformT<T> waveformAbs(n);
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

template bool generateLowResWaveform<TSampleI16>(const TWaveformViewT<TSampleI16> & waveform, TWaveformT<TSampleI16> & waveformLowRes, int nWindow);
