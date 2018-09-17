cmake_minimum_required (VERSION 2.8)
project (KBD-Audio)

set(CMAKE_EXPORT_COMPILE_COMMANDS "on")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Threads REQUIRED)
find_package(SDL2 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(FFTW REQUIRED)

string(STRIP ${SDL2_LIBRARIES} SDL2_LIBRARIES)

if (APPLE)
    find_library(COREFOUNDATION_LIBRARY CoreFoundation )
    find_library(COREGRAPHICS_LIBRARY CoreGraphics )
else (APPLE)
    unset(COREFOUNDATION_LIBRARY)
    unset(COREGRAPHICS_LIBRARY)
endif (APPLE)

include_directories(imgui)
include_directories(imgui/examples)
include_directories(imgui/examples/libs/gl3w)
include_directories(${SDL2_INCLUDE_DIRS})

add_library(Core STATIC
    audio_logger.cpp
    )

target_include_directories(Core PRIVATE
    ${FFTW_INCLUDE_DIRS}
    )

target_link_libraries(Core PRIVATE
    ${CMAKE_THREAD_LIBS_INIT}
    ${SDL2_LIBRARIES}
    ${FFTWF_LIBRARIES}
    ${COREGRAPHICS_LIBRARY}
    ${COREFOUNDATION_LIBRARY}
    )

add_library(ImGui STATIC
    imgui/imgui.cpp
    imgui/imgui_draw.cpp
    imgui/imgui_demo.cpp
    imgui/imgui_widgets.cpp
    imgui/examples/libs/gl3w/GL/gl3w.c
    imgui/examples/imgui_impl_sdl.cpp
    imgui/examples/imgui_impl_opengl3.cpp
    )

target_compile_definitions(ImGui PRIVATE
    IMGUI_USER_CONFIG="../imconfig-vtx32.h"
    )

target_link_libraries(ImGui PRIVATE
    ${CMAKE_DL_LIBS}
    ${SDL2_LIBRARIES}
    ${OPENGL_LIBRARIES}
    ${COREFOUNDATION_LIBRARY}
    )

add_executable(key_average_gui key_average_gui.cpp)
target_link_libraries(key_average_gui PRIVATE Core ImGui)

add_executable(guess_qp guess_qp.cpp)
target_link_libraries(guess_qp PRIVATE Core)

add_executable(guess_qp2 guess_qp2.cpp)
target_link_libraries(guess_qp2 PRIVATE Core)

add_executable(guess_all guess_all.cpp)
target_link_libraries(guess_all PRIVATE Core)

add_executable(key_detector key_detector.cpp)
target_link_libraries(key_detector PRIVATE Core)

add_executable(play play.cpp)
target_link_libraries(play PRIVATE Core)

add_executable(scale scale.cpp)
target_link_libraries(scale PRIVATE Core)

if (APPLE)
    add_executable(record record.cpp key_logger.cpp)
    target_link_libraries(record PRIVATE Core)
endif (APPLE)
