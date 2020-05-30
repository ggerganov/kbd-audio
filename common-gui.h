/*! \file common-gui.h
 *  \brief Enter description here.
 */

#pragma once

#include "common.h"

template<typename T>
bool generateLowResWaveform(const TWaveformViewT<T> & waveform, TWaveformT<T> & waveformLowRes, int nWindow);

template<typename T>
bool generateLowResWaveform(const TWaveformT<T> & waveform, TWaveformT<T> & waveformLowRes, int nWindow) {
    return generateLowResWaveform(getView(waveform, 0), waveformLowRes, nWindow);
}
