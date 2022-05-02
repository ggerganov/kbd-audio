cmake_minimum_required (VERSION 2.8)
project (KBD-Audio)

set(CMAKE_EXPORT_COMPILE_COMMANDS "on")

# options

option(BUILD_EXPERIMENTAL "Build experimental tools" OFF)
option(USE_FINDSDL2 "Use the FindSDL2.cmake script" OFF)

option(KBD_AUDIO_SANITIZE_THREAD         "kbd-audio: enable thread sanitizer" OFF)
option(KBD_AUDIO_SANITIZE_ADDRESS        "kbd-audio: enable address sanitizer" OFF)
option(KBD_AUDIO_SANITIZE_UNDEFINED      "kbd-audio: enable undefined sanitizer" OFF)

# sanitizers

if (KBD_AUDIO_SANITIZE_THREAD)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=thread")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread")
endif()

if (KBD_AUDIO_SANITIZE_ADDRESS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
endif()

if (KBD_AUDIO_SANITIZE_UNDEFINED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=undefined")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined")
endif()

# sdl2

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
if (USE_FINDSDL2)
    set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake/sdl2")
endif()

# build types

include(cmake/GitVars.cmake)
include(cmake/BuildTypes.cmake)

configure_file(${CMAKE_SOURCE_DIR}/build-vars.h.in ${CMAKE_SOURCE_DIR}/build-vars.h @ONLY)

# warnings

if (CMAKE_COMPILER_IS_GNUCC OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
endif()

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

if (EMSCRIPTEN)
    set (CMAKE_CXX_FLAGS " \
    -s USE_SDL=2 \
    -s USE_PTHREADS=1 \
    ")
    #-s DISABLE_EXCEPTION_CATCHING=1 \

    set(CMAKE_EXE_LINKER_FLAGS " \
    --bind \
    --use-preload-cache \
    -s ASSERTIONS=1 \
    -s NO_EXIT_RUNTIME=0 \
    -s PTHREAD_POOL_SIZE=16 \
    -s INITIAL_MEMORY=536870912 \
    ")
elseif(MINGW)
    find_package(PkgConfig REQUIRED)
    pkg_search_module(SDL2 REQUIRED sdl2)

    set(OPENGL_INCLUDE_DIR "")
    set(OPENGL_LIBRARIES "opengl32")
else()
    find_package(OpenGL REQUIRED)
    find_package(SDL2)
    find_package(FFTW)
endif()

find_package(Threads REQUIRED)

if (NOT USE_FINDSDL2 AND NOT SDL2_FOUND AND NOT EMSCRIPTEN)
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
include_directories(${OPENGL_INCLUDE_DIR})

add_library(Core STATIC
    common.cpp
    audio-logger.cpp
    )

target_include_directories(Core PRIVATE
    )

target_link_libraries(Core PRIVATE
    ${CMAKE_THREAD_LIBS_INIT}
    ${SDL2_LIBRARIES}
    ${COREFOUNDATION_LIBRARY}
    )

# todo : this is ugly, what's the proper way?
if (MINGW)
    target_link_libraries(Core PRIVATE stdc++)
endif()

if (EMSCRIPTEN)
    add_library(Gui STATIC
        common-gui.cpp
        imgui/imgui.cpp
        imgui/imgui_draw.cpp
        imgui/imgui_demo.cpp
        imgui/imgui_widgets.cpp
        imgui/examples/imgui_impl_sdl.cpp
        imgui/examples/imgui_impl_opengl3.cpp
        )

    target_compile_definitions(Gui PUBLIC
        IMGUI_USER_CONFIG="../imconfig-vtx32.h"
        )
else()
    add_library(Gui STATIC
        common-gui.cpp
        imgui/imgui.cpp
        imgui/imgui_draw.cpp
        imgui/imgui_demo.cpp
        imgui/imgui_widgets.cpp
        imgui/examples/libs/gl3w/GL/gl3w.c
        imgui/examples/imgui_impl_sdl.cpp
        imgui/examples/imgui_impl_opengl3.cpp
        )

    # force GL3W loader
    target_compile_definitions(Gui PUBLIC
        IMGUI_IMPL_OPENGL_LOADER_GL3W=1
        IMGUI_USER_CONFIG="../imconfig-vtx32.h"
        )

    target_link_libraries(Gui PRIVATE
        ${CMAKE_DL_LIBS}
        ${SDL2_LIBRARIES}
        ${OPENGL_LIBRARIES}
        ${COREFOUNDATION_LIBRARY}
        )
endif()

add_executable(keytap-gui keytap-gui.cpp)
target_link_libraries(keytap-gui PRIVATE Core Gui)

add_executable(keytap2-gui keytap2-gui.cpp subbreak2.cpp)
target_link_libraries(keytap2-gui PRIVATE Core Gui)

add_executable(keytap3-gui keytap3-gui.cpp subbreak3.cpp)
target_link_libraries(keytap3-gui PRIVATE Core Gui)

add_executable(keytap3-app keytap3-app.cpp subbreak3.cpp)
target_link_libraries(keytap3-app PRIVATE Core)

add_executable(view-full-gui view-full-gui.cpp)
target_link_libraries(view-full-gui PRIVATE Core Gui)

if (EMSCRIPTEN)
    # view-full-gui
    set(TARGET view-full-gui)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${TARGET})

    set_target_properties(${TARGET} PROPERTIES LINK_FLAGS " \
    -s TOTAL_MEMORY=536870912 \
    -s LZ4=1 \
    --preload-file ${PROJECT_SOURCE_DIR}/chal/record@/chal/record/ \  \
    ")

    ## keytap-gui
    set(TARGET keytap-gui)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${TARGET})

    set_target_properties(${TARGET} PROPERTIES LINK_FLAGS " \
    -s TOTAL_MEMORY=536870912 \
    -s FORCE_FILESYSTEM=1 \
    ")

    # keytap2-gui
    set(TARGET keytap2-gui)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${TARGET})

    set_target_properties(${TARGET} PROPERTIES LINK_FLAGS " \
    -s TOTAL_MEMORY=536870912 \
    -s FORCE_FILESYSTEM=1 \
    -s LZ4=1 \
    --preload-file ${PROJECT_SOURCE_DIR}/data/@/data/ \
    ")

    configure_file(${CMAKE_SOURCE_DIR}/index-${TARGET}-tmpl.html ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/index.html @ONLY)
    configure_file(${CMAKE_SOURCE_DIR}/style.css ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/style.css @ONLY)

    # keytap3-app
    set(TARGET keytap3-app)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${TARGET})

    set_target_properties(${TARGET} PROPERTIES LINK_FLAGS " \
    -s FORCE_FILESYSTEM=1 \
    -s LZ4=1 \
    --preload-file ${PROJECT_SOURCE_DIR}/data/ggwords-6-gram.dat.binary@/data/ \
    ")

    configure_file(${CMAKE_SOURCE_DIR}/index-${TARGET}-tmpl.html ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/index.html @ONLY)
    configure_file(${CMAKE_SOURCE_DIR}/style.css ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/style.css @ONLY)

    # keytap3-gui
    set(TARGET keytap3-gui)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${TARGET})

    set_target_properties(${TARGET} PROPERTIES LINK_FLAGS " \
    -s TOTAL_MEMORY=536870912 \
    -s FORCE_FILESYSTEM=1 \
    -s LZ4=1 \
    --preload-file ${PROJECT_SOURCE_DIR}/data/ggwords-6-gram.dat.binary@/data/ \
    ")

    configure_file(${CMAKE_SOURCE_DIR}/index-${TARGET}-tmpl.html ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/index.html @ONLY)
    configure_file(${CMAKE_SOURCE_DIR}/style.css ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/style.css @ONLY)
endif()

if (NOT EMSCRIPTEN)
    add_executable(key-detector key-detector.cpp)
    target_link_libraries(key-detector PRIVATE Core)

    add_executable(keytap keytap.cpp)
    target_link_libraries(keytap PRIVATE Core)

    add_executable(keytap2 keytap2.cpp)
    target_link_libraries(keytap2 PRIVATE Core)

    add_executable(keytap3 keytap3.cpp subbreak3.cpp)
    target_link_libraries(keytap3 PRIVATE Core)

    add_executable(play play.cpp)
    target_link_libraries(play PRIVATE Core)

    add_executable(play-full play-full.cpp)
    target_link_libraries(play-full PRIVATE Core)

    if (NOT MINGW)
        # todo : need to find termios.h alternative for Windows
        add_executable(record record.cpp)
        target_link_libraries(record PRIVATE Core)
    endif()

    add_executable(record-full record-full.cpp)
    target_link_libraries(record-full PRIVATE Core)

    add_executable(view-gui view-gui.cpp)
    target_link_libraries(view-gui PRIVATE Core Gui)

    add_executable(compress-n-grams compress-n-grams.cpp subbreak3.cpp)
    target_link_libraries(compress-n-grams PRIVATE Core)

    #
    ## Experimental stuff

    if (BUILD_EXPERIMENTAL)
        if (FFTW_FOUND)
            add_executable(key-average-gui key-average-gui.cpp)
            target_include_directories(key-average-gui PRIVATE ${FFTW_INCLUDE_DIRS})
            target_link_libraries(key-average-gui PRIVATE Core Gui ${FFTW_LIBRARIES})
        else()
            message(WARNING "Skipping 'key-average-gui' target because FFTW is not available")
        endif()

        add_executable(keytap3-multi keytap3-multi.cpp subbreak3.cpp)
        target_link_libraries(keytap3-multi PRIVATE Core)

        add_executable(guess-qp guess-qp.cpp)
        target_link_libraries(guess-qp PRIVATE Core)

        add_executable(guess-qp2 guess-qp2.cpp)
        target_link_libraries(guess-qp2 PRIVATE Core)

        add_executable(scale scale.cpp)
        target_link_libraries(scale PRIVATE Core)

        add_executable(subbreak subbreak.cpp)
        target_link_libraries(subbreak PRIVATE Core)

        add_executable(non-exact-subbreak non-exact-subbreak.cpp subbreak2.cpp)
        target_link_libraries(non-exact-subbreak PRIVATE Core)

        add_executable(non-exact-subbreak2 non-exact-subbreak2.cpp subbreak2.cpp)
        target_link_libraries(non-exact-subbreak2 PRIVATE Core)

        add_executable(generate-clusters generate-clusters.cpp subbreak2.cpp)
        target_link_libraries(generate-clusters PRIVATE Core)

        add_executable(test-subbreak3 test-subbreak3.cpp subbreak3.cpp)
        target_link_libraries(test-subbreak3 PRIVATE Core)
    endif()

endif()
