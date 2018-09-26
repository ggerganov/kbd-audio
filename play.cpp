/*! \file play.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include <SDL2/SDL.h>
#include <SDL2/SDL_audio.h>

#include "audio_logger.h"

#include <fstream>

bool g_terminate = false;

void cbPlayback(void * userdata, uint8_t * stream, int len) {
    std::ifstream * fin = (std::ifstream *)(userdata);
    int keyPressed = -1;
    fin->read((char *)(&keyPressed), sizeof(keyPressed));
    if (fin->eof()) {
        g_terminate = true;
        return;
    }
    printf("%c", keyPressed);
    fflush(stdout);
    fin->read((char *)(stream), len); // todo
    if (fin->eof()) g_terminate = true;
}

int main(int argc, char ** argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s input.kbd\n", argv[0]);
        return -127;
    }

    std::ifstream fin(argv[1], std::ios::binary);

    int bufferSize_frames = 1;
    fin.read((char *)(&bufferSize_frames), sizeof(bufferSize_frames));
    printf("Buffer size = %d frames\n", bufferSize_frames);

    if (SDL_Init(SDL_INIT_AUDIO) < 0) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't initialize SDL: %s\n", SDL_GetError());
        return -1;
    }

    int nDevices = SDL_GetNumAudioDevices(SDL_FALSE);
    printf("Found %d playback devices:\n", nDevices);
    for (int i = 0; i < nDevices; i++) {
        printf("    - Playback device #%d: '%s'\n", i, SDL_GetAudioDeviceName(i, SDL_FALSE));
    }

    SDL_AudioSpec playbackSpec;
    SDL_zero(playbackSpec);

    playbackSpec.freq = 24000;
    playbackSpec.format = AUDIO_F32SYS;
    playbackSpec.channels = 1;
    playbackSpec.samples = bufferSize_frames*AudioLogger::kSamplesPerFrame;
    playbackSpec.callback = cbPlayback;
    playbackSpec.userdata = (void *)(&fin);

    SDL_AudioSpec obtainedSpec;
    SDL_zero(obtainedSpec);

    auto deviceIdOut = SDL_OpenAudioDevice(NULL, SDL_FALSE, &playbackSpec, &obtainedSpec, 0);
    if (!deviceIdOut) {
        SDL_LogError(SDL_LOG_CATEGORY_APPLICATION, "Couldn't open an audio device for playback: %s!\n", SDL_GetError());
        SDL_Quit();
        return -2;
    }

    int sampleSize_bytes = 4; // todo

    printf("Opened playback device %d\n", deviceIdOut);
    printf("    Frequency:  %d\n", obtainedSpec.freq);
    printf("    Format:     %d (%d bytes)\n", obtainedSpec.format, sampleSize_bytes);
    printf("    Channels:   %d\n", obtainedSpec.channels);
    printf("    Samples:    %d\n", obtainedSpec.samples);

    SDL_PauseAudioDevice(deviceIdOut, 0);

    while(g_terminate == false) {
        SDL_Delay(100);
    }

    fin.close();

    SDL_CloseAudio();

    return 0;
}
