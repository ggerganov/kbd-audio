cmake_minimum_required (VERSION 2.8)
project (KBD-Audio)

option(BUILD_EXPERIMENTAL "Build experimental tools" OFF)
option(USE_FINDSDL2 "Use the FindSDL2.cmake script" OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS "on")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
if (USE_FINDSDL2)
    set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake/sdl2")
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "DEBUG" "RELEASE" "RELWITHDEBINFO")
endif()

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Threads REQUIRED)
find_package(OpenGL REQUIRED)
find_package(SDL2)
find_package(FFTW)

if (NOT USE_FINDSDL2 AND NOT SDL2_FOUND)
    message(WARNING "Unable to find SDL2 library. It is either not installed or CMake cannot find it."
        " In the latter case, setting the USE_FINDSDL2 variable might help:\n"
        "   $ cmake -D USE_FINDSDL2 .."
        )

    message(FATAL_ERROR "Aborting")
endif()

if (NOT FFTW_FOUND)
    message(WARNING "FFTW library not available. Some targets will not be built")
endif()

string(STRIP "${SDL2_LIBRARIES}" SDL2_LIBRARIES)

if (APPLE)
    find_library(COREGRAPHICS_LIBRARY CoreGraphics )
    find_library(COREFOUNDATION_LIBRARY CoreFoundation )
else (APPLE)
    unset(COREGRAPHICS_LIBRARY)
    unset(COREFOUNDATION_LIBRARY)
endif (APPLE)

include_directories(imgui)
include_directories(imgui/examples)
include_directories(imgui/examples/libs/gl3w)
include_directories(${SDL2_INCLUDE_DIRS})

add_library(Core STATIC
    audio_logger.cpp
    )

target_include_directories(Core PRIVATE
    )

target_link_libraries(Core PRIVATE
    ${CMAKE_THREAD_LIBS_INIT}
    ${SDL2_LIBRARIES}
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

add_executable(keytap keytap.cpp)
target_link_libraries(keytap PRIVATE Core)

add_executable(keytap2 keytap2.cpp)
target_link_libraries(keytap2 PRIVATE Core)

add_executable(keytap-gui keytap-gui.cpp)
target_link_libraries(keytap-gui PRIVATE Core ImGui)

add_executable(keytap2-gui keytap2-gui.cpp subbreak2.cpp)
target_link_libraries(keytap2-gui PRIVATE Core ImGui)

add_executable(play play.cpp)
target_link_libraries(play PRIVATE Core)

add_executable(play-full play-full.cpp)
target_link_libraries(play-full PRIVATE Core)

add_executable(record record.cpp)
target_link_libraries(record PRIVATE Core)

add_executable(record-full record-full.cpp)
target_link_libraries(record-full PRIVATE Core)

add_executable(view-gui view-gui.cpp)
target_link_libraries(view-gui PRIVATE Core ImGui)

add_executable(view-full-gui view-full-gui.cpp)
target_link_libraries(view-full-gui PRIVATE Core ImGui)

#
## Experimental stuff

if (BUILD_EXPERIMENTAL)
    if (FFTW_FOUND)
        add_executable(key_average_gui key_average_gui.cpp)
        target_include_directories(key_average_gui PRIVATE ${FFTW_INCLUDE_DIRS})
        target_link_libraries(key_average_gui PRIVATE Core ImGui ${FFTW_LIBRARIES})
    else()
        message(WARNING "Skipping 'key_average_gui' target because FFTW is not available")
    endif()

    add_executable(guess_qp guess_qp.cpp)
    target_link_libraries(guess_qp PRIVATE Core)

    add_executable(guess_qp2 guess_qp2.cpp)
    target_link_libraries(guess_qp2 PRIVATE Core)

    add_executable(key_detector key_detector.cpp)
    target_link_libraries(key_detector PRIVATE Core)

    add_executable(scale scale.cpp)
    target_link_libraries(scale PRIVATE Core)

    add_executable(subbreak subbreak.cpp)
    target_link_libraries(subbreak PRIVATE Core)
endif()
