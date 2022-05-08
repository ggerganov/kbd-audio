kbd-audio
=====
[![Actions Status](https://github.com/ggerganov/kbd-audio/workflows/CI/badge.svg)](https://github.com/ggerganov/kbd-audio/actions)

This is a collection of command-line and GUI tools for capturing and analyzing audio data.

### Keytap

The most interesting tool is called **keytap** - it can guess pressed keyboard keys only by analyzing the audio captured from the computer's microphone.

Check this blog post for more details:

[Keytap: description and some random thoughts](https://ggerganov.github.io/jekyll/update/2018/11/30/keytap-description-and-thoughts.html)

[Video: short demo of Keytap in action](https://www.youtube.com/watch?v=2OjzI9m7W10)

[Try it online:](https://ggerganov.github.io/keytap)

<a href="https://ggerganov.github.io/keytap" target="_blank"><img src="https://i.imgur.com/FXa60Pr.gif" style="display: inline-block; overflow: hidden; width: 99%;"></img></a>

### Keytap2

The **keytap2** tool is another interesting tool for recovering text from audio. It does not require training data - instead it uses statistical information about the frequencies of the letters and n-grams in the English language.

A more detailed description of the tool is available here: [Keytap2 discussion](https://github.com/ggerganov/kbd-audio/discussions/31)

[Video: short demo of Keytap2 in action](https://www.youtube.com/watch?v=jNtw17S6SR0)

[CTF: can you guess the text being typed?](https://ggerganov.github.io/keytap-challenge/)

[Try it online:](https://keytap2.ggerganov.com)

<a href="https://keytap2.ggerganov.com" target="_blank"><img src="https://i.imgur.com/nPlLEDN.jpg" style="display: inline-block; overflow: hidden; width: 99%;"></img></a>

### Keytap3

This version introduces significant algorithm improvements and better n-gram statistics compared to keytap2. The attack is now fully
automated and does not require any manual intervation during the text recovery process.

[Video: short demo of using Keytap3](https://youtu.be/5aphvxpSt3o)

[Video: another example of using Keytap3](https://youtu.be/kCOrxrR-4ak)

[GUI for Keytap3](https://keytap3-gui.ggerganov.com)

[Check if your keyboard is vulnerable to Keytap:](https://keytap3.ggerganov.com)

<a href="https://keytap3.ggerganov.com" target="_blank"><img src="https://user-images.githubusercontent.com/1991296/166096331-ab26f7f8-08e0-48d6-abd7-57017ebf1866.JPEG" style="display: inline-block; overflow: hidden; width: 99%;"></img></a>

### What people say about Keytap

*"This works incredibly well.\
I hope you realize what you've created (and made available to every person in the world)."* -- ffpip

*"I just tried it and it works incredibly well. It kind of makes me want to stop using a mechanical keyboard."* -- Karawebnetwork

*"This attack and Van Eck phreaking are why Edward Snowden, while typing passwords and other sensitive information, would pull a blanket over himself and his laptop."* -- aarchi

*"This is what mechanical keyboard users deserve"* -- super guy

*"fuck.."* -- Lluis Franco

## Build instructions

Dependencies:

 - **SDL2** - used to capture audio and to open GUI windows [libsdl](https://www.libsdl.org)

       [Ubuntu]
       $ sudo apt install libsdl2-dev

       [Mac OS with brew]
       $ brew install sdl2

       [MSYS2]
       $ pacman -S git cmake make mingw-w64-x86_64-dlfcn mingw-w64-x86_64-gcc mingw-w64-x86_64-SDL2

 - **FFTW3** *(optional)* - some of the helper tools perform Fourier transformations [fftw](http://www.fftw.org)

**Linux, FreeBSD, Mac OS, Windows (MSYS2 + MinGW)**

    git clone https://github.com/ggerganov/kbd-audio
    cd kbd-audio
    git submodule update --init
    mkdir build && cd build
    cmake ..
    make

## Tools

Short summary of the available tools. If the status of the tool is not **stable**, expect problems and non-optimal results.

| Name                | Type    | Status      |
| ---                 | ---     | ---         |
| **record**          | text    | **stable**  |
| **record-full**     | text    | **stable**  |
| **play**            | text    | **stable**  |
| **play-full**       | text    | **stable**  |
| **view-gui**        | gui     | **stable**  |
| **view-full-gui**   | gui     | **stable**  |
| **key-detector**    | text    | **stable**  |
| **keytap**          | text    | **stable**  |
| **keytap-gui**      | gui     | **stable**  |
| **keytap2-gui**     | gui     | **stable**  |
| **keytap3**         | text    | **stable**  |
| **keytap3-gui**     | gui     | **stable**  |
| -                   | *extra* | -           |
| **guess-qp**        | text    | experiment  |
| **guess-qp2**       | text    | experiment  |
| **keytap3-multi**   | text    | experiment  |
| **scale**           | text    | experiment  |
| **subreak**         | text    | experiment  |
| **key-average-gui** | gui     | experiment  |
| **keytap2**         | text    | experiment  |

## Tool details

* **record-full**

  Record audio to a raw binary file on disk

      ./record-full output.kbd [-cN]

  ---

* **play-full**

  Playback a recording captured via the **record-full** tool

      ./play-full input.kbd [-pN]

  ---

* **record**

  Record audio only while typing. Useful for collecting training data for **keytap**

      ./record output.kbd [-cN] [-CN]

  ---

* **play**

  Playback a recording created via the **record** tool

      ./play input.kbd [-pN]

  ---

* **keytap**

  Detect pressed keys via microphone audio capture in real-time. Uses training data captured via the **record** tool.

      ./keytap input0.kbd [input1.kbd] [input2.kbd] ... [-cN] [-CN] [-pF] [-tF]

  ---

* **keytap-gui**

  Detect pressed keys via microphone audio capture in real-time. Uses training data captured via the **record** tool. GUI version.

      ./keytap-gui input0.kbd [input1.kbd] [input2.kbd] ... [-cN] [-CN]

  Online demo: https://keytap.ggerganov.com

  ---

* **keytap2-gui** record.kbd n-gram-dir [-pN] [-cN] [-CN]

  Detect pressed keys via microphone audio capture. Uses statistical information (n-gram frequencies) about the language. **No training data is required**. The *'record.kbd'* input file has to be generated via the **record-full** tool and contains the audio data that will be analyzed. The *'n-gram-dir'* folder file has to contain n-gram probability files for the corresponding language.

      ./keytap2-gui record.kbd ../data

  Online demo: https://keytap2.ggerganov.com

  ---

* **keytap3**

  Fully automated recovery of unknown text from audio recordings.

      ./keytap3 input.kbd ../data [-cN] [-CN] [-pF] [-tF] [-FN] [-fN]

  Online demo: https://keytap3.ggerganov.com

  ---

* **keytap3-gui**

  GUI version of the **keytap3** tool.

      ./keytap3-gui input.kbd ../data [-cN] [-CN] [-pF] [-tF] [-FN] [-fN]

  Online demo: https://keytap3-gui.ggerganov.com

  ---

* **view-full-gui**

  Visualize waveforms recorded with the **record-full** tool. Can also playback the audio data.

      ./view-full-gui input.kbd [-pN]

  <a href="https://i.imgur.com/scjTaXw.png" target="_blank">![view-full-gui](https://i.imgur.com/scjTaXw.png)</a>

  ---

* **view-gui**

  Visualize training data recorded with the **record** tool. Can also playback the audio data.

      ./view-gui input.kbd [-pN]

  <a href="https://i.imgur.com/2binGaZ.png" target="_blank">![view-full-gui](https://i.imgur.com/2binGaZ.png)</a>

  ---

## Feedback

Any feedback about the performance of the tools is highly appreciated. Please drop a comment [here](https://github.com/ggerganov/kbd-audio/issues/3).
