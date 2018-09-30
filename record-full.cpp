/*! \file record-full.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "constants.h"
#include "audio_logger.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_audio.h>

#include <fstream>

bool g_terminate = false;

void cbPlayback(void * userdata, uint8_t * stream, int len) {
    std::ofstream * fout = (std::ofstream *)(userdata);
    fout->write((char *)(stream), len); // todo
}

int main(int argc, const char ** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s output.kbd\n", argv[0]);
        return -127;
    }

    std::ofstream fout(argv[1], std::ios::binary);
    if (fout.good() == false) {
        fprintf(stderr, "Failed to open file '%s'\n", argv[1]);
        return -1;
    }

    int bufferSize_frames = 1;

    if (SDL_Init(SDL_INIT_AUDIO) < 0) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't initialize SDL: %s\n", SDL_GetError());
        return -1;
    }

    int nDevices = SDL_GetNumAudioDevices(SDL_TRUE);
    printf("Found %d capture devices:\n", nDevices);
    for (int i = 0; i < nDevices; i++) {
        printf("    - Capture device #%d: '%s'\n", i, SDL_GetAudioDeviceName(i, SDL_TRUE));
    }

    SDL_AudioSpec captureSpec;
    SDL_zero(captureSpec);

    captureSpec.freq = 24000;
    captureSpec.format = AUDIO_F32SYS;
    captureSpec.channels = 1;
    captureSpec.samples = AudioLogger::kSamplesPerFrame;
    captureSpec.callback = cbPlayback;
    captureSpec.userdata = (void *)(&fout);

    SDL_AudioSpec obtainedSpec;
    SDL_zero(obtainedSpec);

    auto deviceIdIn = SDL_OpenAudioDevice(NULL, SDL_TRUE, &captureSpec, &obtainedSpec, 0);
    if (!deviceIdIn) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't open an audio device for capture: %s!\n", SDL_GetError());
        SDL_Quit();
        return -2;
    }

    int sampleSize_bytes = 4; // todo

    printf("Opened capture device %d\n", deviceIdIn);
    printf("    Frequency:  %d\n", obtainedSpec.freq);
    printf("    Format:     %d (%d bytes)\n", obtainedSpec.format, sampleSize_bytes);
    printf("    Channels:   %d\n", obtainedSpec.channels);
    printf("    Samples:    %d\n", obtainedSpec.samples);

    SDL_PauseAudioDevice(deviceIdIn, 0);

    while(g_terminate == false) {
        SDL_Delay(100);
    }

    fout.close();

    SDL_CloseAudio();

    return 0;
}
