#!/bin/bash

    #-s WASM=1 \
    #-s ALLOW_MEMORY_GROWTH=1 \
    #-s USE_SDL=2 \

    #-s USE_PTHREADS=1 \
    #-s PTHREAD_POOL_SIZE=2 \

    #./main-sdl.cpp \

    #../../imgui/examples/imgui_impl_sdl.cpp \
    #../../imgui/examples/imgui_impl_glfw.cpp \

em++ -O3 -std=c++17 \
    -DCG_OUTPUT -DCG_LOGGING=10 \
    -s WASM=1 \
    -s ASSERTIONS=1 \
    -s GL_ASSERTIONS=1 \
    -s GL_DEBUG=1 \
    -s TOTAL_MEMORY=536870912 \
    -s USE_PTHREADS=1 \
    -s PTHREAD_POOL_SIZE=8 \
    -s USE_SDL=2 \
    -s EXPORTED_FUNCTIONS='["_keyPressedCallback", "_main"]' \
    -s EXTRA_EXPORTED_RUNTIME_METHODS='["ccall", "cwrap"]' \
    -I ../.. \
    -I ../../imgui \
    -I ../../imgui/examples \
    -I ../../imgui/examples/libs/gl3w \
    ../../common.cpp \
    ../../common-gui.cpp \
    ../../keytap-gui.cpp \
    ../../audio_logger.cpp \
    ../../imgui/imgui.cpp \
    ../../imgui/imgui_draw.cpp \
    ../../imgui/imgui_demo.cpp \
    ../../imgui/imgui_widgets.cpp \
    ../../imgui/examples/imgui_impl_sdl.cpp \
    ../../imgui/examples/imgui_impl_opengl3.cpp \
    -o keytap-gui.js

cp -v keytap-gui.js         ../../../ggerganov.github.io/scripts/keytap-gui/
cp -v keytap-gui.wasm       ../../../ggerganov.github.io/scripts/keytap-gui/
cp -v keytap-gui.worker.js  ../../../ggerganov.github.io/scripts/keytap-gui/
