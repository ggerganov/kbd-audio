/*! \file record-full.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include "constants.h"
#include "common.h"

#include <SDL.h>
#include <SDL_audio.h>

#include <fstream>

bool g_terminate = false;

void cbPlayback(void * userdata, uint8_t * stream, int len) {
    std::ofstream * fout = (std::ofstream *)(userdata);
    fout->write((char *)(stream), len); // todo
}

int main(int argc, char ** argv) {
    printf("Usage: %s output.kbd [-cN]\n", argv[0]);
    printf("    -cN - select capture device N\n");
    printf("\n");

    if (argc < 2) {
        return -127;
    }

    auto argm = parseCmdArguments(argc, argv);
    int captureId = argm["c"].empty() ? 0 : std::stoi(argm["c"]);

    std::ofstream fout(argv[1], std::ios::binary);
    if (fout.good() == false) {
        fprintf(stderr, "Failed to open file '%s'\n", argv[1]);
        return -1;
    }

    if (SDL_Init(SDL_INIT_AUDIO) < 0) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't initialize SDL: %s\n", SDL_GetError());
        return -1;
    }

    int nDevices = SDL_GetNumAudioDevices(SDL_TRUE);
    printf("Found %d capture devices:\n", nDevices);
    for (int i = 0; i < nDevices; i++) {
        printf("    - Capture device #%d: '%s'\n", i, SDL_GetAudioDeviceName(i, SDL_TRUE));
    }

    if (captureId < 0 || captureId >= nDevices) {
        printf("Invalid capture device id selected - %d\n", captureId);
        return -1;
    }

    SDL_AudioSpec captureSpec;
    SDL_zero(captureSpec);

    captureSpec.freq = kSampleRate;
    captureSpec.format = AUDIO_F32SYS;
    captureSpec.channels = 1;
    captureSpec.samples = kSamplesPerFrame;
    captureSpec.callback = cbPlayback;
    captureSpec.userdata = (void *)(&fout);

    SDL_AudioSpec obtainedSpec;
    SDL_zero(obtainedSpec);

    printf("Attempt to open capture device %d : '%s' ...\n", captureId, SDL_GetAudioDeviceName(captureId, SDL_TRUE));
    auto deviceIdIn = SDL_OpenAudioDevice(SDL_GetAudioDeviceName(captureId, SDL_TRUE), SDL_TRUE, &captureSpec, &obtainedSpec, 0);
    if (!deviceIdIn) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't open an audio device for capture: %s!\n", SDL_GetError());
        SDL_Quit();
        return -2;
    }

    int sampleSize_bytes = 4; // todo

    printf("Opened capture device succesfully!\n");
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
