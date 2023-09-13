# BASS Ginsu Synthesizer Library

This is a Ginsu synthesizer made mostly with BASS.

A small portion uses vgmstream to load in the sample from a Ginsu file.



This is not designed to directly play back Ginsu files.

This is designed to give a signed 16-bit PCM stream with the GetData function of a Ginsu stream which you would plug to your own sound engine.

## Basic usage

1. Set up a project with BASS, BASSmix and vgmstream (make sure their .lib files are included in the linker)

2. Include "BASSGinsu.hpp"

3. Read the top of "BASSGinsu.hpp" for further instructions

## TODO

There is still a lot to be done. Check the top of the source file for more information.


