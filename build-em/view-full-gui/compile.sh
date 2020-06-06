#!/bin/bash

em++ -O3 -std=c++17 \
    -DCG_OUTPUT -DCG_LOGGING=10 \
    -s WASM=1 \
    -s ASSERTIONS=1 \
    -s TOTAL_MEMORY=134217728 \
    -s USE_SDL=2 \
    -s EXPORTED_FUNCTIONS='["_doInit", "_main"]' \
    -s EXTRA_EXPORTED_RUNTIME_METHODS='["ccall", "cwrap"]' \
    --preload-file ../../chal/record \
    -I ../.. \
    -I ../../imgui \
    -I ../../imgui/examples \
    -I ../../imgui/examples/libs/gl3w \
    ../../common.cpp \
    ../../common-gui.cpp \
    ../../view-full-gui.cpp \
    ../../imgui/imgui.cpp \
    ../../imgui/imgui_draw.cpp \
    ../../imgui/imgui_demo.cpp \
    ../../imgui/imgui_widgets.cpp \
    ../../imgui/examples/imgui_impl_sdl.cpp \
    ../../imgui/examples/imgui_impl_opengl3.cpp \
    -o view-full-gui.js

cp -v view-full-gui.js     ../../../ggerganov.github.io/scripts/view-full-gui/
#cp -v view-full-gui.js.mem ../../../ggerganov.github.io/scripts/view-full-gui/
cp -v view-full-gui.data ../../../ggerganov.github.io/scripts/view-full-gui/
cp -v view-full-gui.wasm   ../../../ggerganov.github.io/scripts/view-full-gui/
#cp -v pthread-main.js   ../../../ggerganov.github.io/scripts/view-full-gui/
