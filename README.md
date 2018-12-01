# kbd-audio [![Build Status](https://travis-ci.org/ggerganov/kbd-audio.svg?branch=master)](https://travis-ci.org/ggerganov/kbd-audio?branch=master)

## Description

This is a collection of command-line and GUI tools for capturing and analyzing audio data. The most interesting tool is called **keytap** - it can guess pressed keyboard keys only by analyzing the audio captured from the computer's microphone.

Check this blog post for more details:

[Keytap: description and some random thoughts](https://ggerganov.github.io/jekyll/update/2018/11/30/keytap-description-and-thoughts.html)

## Build instructions

Dependencies:

 - **SDL2** - used to capture audio and to open GUI windows [libsdl](https://www.libsdl.org)
 - **FFTW3** - some of the helper tools perform Fourier transformations [fftw](http://www.fftw.org)
 
**Linux and Mac OS**

    git clone https://github.com/ggerganov/kbd-audio
    cd kbd-audio
    git submodule update --init
    mkdir build && cd build
    cmake ..
    make
    
**Windows**

    (todo, PRs welcome)

## Tools

### record-full

Record audio to a raw binary file on disk

**Usage:** ./record-full output.kbd


### play-full

Playback a recording captured via the **record-full** tool

**Usage:** ./play-full input.kbd


### record

Record audio only while typing. Useful for collecting training data for **keytap**

**Usage:** ./record output.kbd


### play

Playback a recording created via the **record** tool

**Usage:** ./play input.kbd


### keytap

Detect pressed keys via microphone audio capture in real-time. Uses training data captured via the **record** tool.

**Usage:** ./keytap-gui input0.kbd [input1.kbd] [input2.kbd] ...

[**Live demo *(WebAssembly threads required)***](https://ggerganov.github.io/jekyll/update/2018/11/24/keytap.html)

<a href="https://i.imgur.com/mnRvT1X.gif" target="_blank">![Keytap](https://i.imgur.com/FXa60Pr.gif)</a>


### keytap2 *(work in progress)*

Detect pressed keys via microphone audio capture. Uses statistical information (n-gram frequencies) about the language. **No training data is required**. The *'recording.kbd'* input file has to be generated via the **record-full** tool and contains the audio data that will be analyzed. The *'n-gram.txt'* file has to contain n-gram probabilities for the corresponding language. 

**Usage:** ./keytap2-gui recording.kbd n-gram.txt

<a href="https://i.imgur.com/yR3m5Bm.jpg" target="_blank">![Keytap](https://i.imgur.com/yR3m5Bm.jpg)</a>

## Feedback

Any feedback about the performance of the tools is highly appreciated. Please drop a comment [here](https://github.com/ggerganov/kbd-audio/issues/3).
