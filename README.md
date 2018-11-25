## record-full

Record audio to a raw binary file on disk

**Usage:** ./record-full output.kbd


## play-full

Playback a recording captured via the **record-full** tool

**Usage:** ./play-full input.kbd


## record

Record audio only while typing. Useful for collecting training data for **keytap**

**Usage:** ./record output.kbd


## play

Playback a recording created via the **record** tool

**Usage:** ./play input.kbd


## keytap

Detect pressed keys via microphone audio capture in real-time. Uses training data captured via the **record** tool.

**Usage:** ./keytap-gui input0.kbd [input1.kbd] [input2.kbd] ...

[**Live demo *(WebAssembly threads required)***](https://ggerganov.github.io/jekyll/update/2018/11/17/wave-gui.html)

<a href="https://i.imgur.com/mnRvT1X.gif" target="_blank">![Keytap](https://i.imgur.com/FXa60Pr.gif)</a>


## keytap2 *(work in progress)*

Detect pressed keys via microphone audio capture. Uses statistical information (n-gram frequencies) about the English language. No training data required. The *'recording.kbd'* input file has to be generated via the **record-full** tool and contains the audio data that will be analyzed.

**Usage:** ./keytap2-gui recording.kbd n-gram.txt

<a-href="https://i.imgur.com/yR3m5Bm.jpg" target="_blank">![Keytap](https://i.imgur.com/yR3m5Bm.jpg)</a>
