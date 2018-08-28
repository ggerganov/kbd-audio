#!/bin/bash

g++ -std=c++14 -o record record.cpp key_logger.cpp audio_logger.cpp -framework CoreFoundation -framework CoreGraphics -framework SDL2
g++ -std=c++14 -o play play.cpp -framework CoreFoundation -framework CoreGraphics -framework SDL2
