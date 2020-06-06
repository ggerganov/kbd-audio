/*! \file common-gui.h
 *  \brief Enter description here.
 */

#pragma once

#include "common.h"

#include <SDL.h>

namespace Gui {
struct Objects {
    SDL_Window * window = nullptr;
    SDL_GLContext gl_context = nullptr;
};

bool init(const char * windowTitle, int windowSizeX, int windowSizeY, Objects & objects);
bool render(const Objects & objects);
bool free(const Objects & objects);
}
