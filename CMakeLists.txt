cmake_minimum_required (VERSION 2.8)
project (KBD-Audio)

set(CMAKE_EXPORT_COMPILE_COMMANDS "on")

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Threads REQUIRED)
find_package(SDL2 REQUIRED)
find_package(OpenGL REQUIRED)

if (APPLE)
    find_library(COREFOUNDATION_LIBRARY CoreFoundation )
    find_library(COREGRAPHICS_LIBRARY CoreGraphics )
endif (APPLE)

include_directories(imgui)
include_directories(imgui/examples)
include_directories(imgui/examples/libs/gl3w)
include_directories(${SDL2_INCLUDE_DIRS})

add_library(Core STATIC
    key_logger.cpp
    audio_logger.cpp
    )

target_link_libraries(Core PRIVATE
    ${CMAKE_THREAD_LIBS_INIT}
    ${SDL2_LIBRARIES}
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

target_link_libraries(ImGui PRIVATE
    ${SDL2_LIBRARIES}
    ${OPENGL_LIBRARIES}
    ${COREFOUNDATION_LIBRARY}
    )

add_executable(key_average_gui key_average_gui.cpp)
target_link_libraries(key_average_gui PRIVATE Core ImGui)

add_executable(record record.cpp)
target_link_libraries(record PRIVATE Core)

add_executable(play play.cpp)
target_link_libraries(play PRIVATE Core)
