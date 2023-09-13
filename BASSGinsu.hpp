//
// BASS Ginsu Synthesizer Library
// by Tenjoin / xan1242
// 
// Requirements:
// - BASS
// - BASSmix
// - vgmstream
// 
// Usage:
// 1. Initialize the player (BASSGinsuPlayer) (creation args are the same as BASS_Init)
// 2. Create a stream with BASSGinsuPlayer::CreateStream for each car. Use a Ginsu file directly.
// 3. Get the signed 16-bit PCM data with BASSGinsu(Multi)Stream::GetData and pipe it to the game sound engine
// 4. Control the streams for each car with BASSGinsu(Multi)Stream's public functions (e.g. to set the RPM use SetFrequency)
// 5. Once done, destroy the BASSGinsuPlayer (or - for finer grain control, release each stream with BASSGinsuPlayer::ReleaseStream for each released car)
//

//
// TODO list:
// - Implement exceptions and error handling!
// - Implement reading of a Ginsu from memory!
// - Tri-crossfading
// - Channel generation on the fly (instead of keeping them all active at once)
// - Redline & idle re-pitch -- done? Maybe this was referring to re-pitching based on the RPM...
// - LPF based on throttle -- maybe implement BASSfx for this? or use the dx8 effect? -- done? not quite LPF but using a param EQ
// - Decel amount based on throttle -- basically, the xfade will depend on how hard the throttle is pushed -- done by tracking the rate?
// - Maybe some control for clutch release and LPF for it too
// - DOCUMENTATION IN GENERAL!
// 
// - ...and other TODOs scattered around the code...
//

#pragma once
#ifndef MOD_BASSGINSU_H
#define MOD_BASSGINSU_H

#define NOMINMAX
#include <iostream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include "FrameLimiter.hpp"
#include <bass.h>
#include <bassmix.h>

#ifdef __cplusplus
extern "C"
{
#endif
#include <vgmstream.h>
#ifdef __cplusplus
}
#endif

#ifdef _MSC_VER
#pragma comment(lib, "bass.lib")
#pragma comment(lib, "bassmix.lib")
#pragma comment(lib, "libvgmstream.lib")
#endif

struct GinsuDataLayout
{
    uint8_t id[4];
    uint8_t ver[2];
    uint16_t flags;
    float minFrequency;
    float maxFrequency;
    uint32_t segCount;
    uint32_t cycleCount;
    uint32_t sampleCount;
    uint32_t sampleRate;
};

enum GIN_DIRECTION
{
    GIN_BACKWARDS = -1,
    GIN_NOWHERE,
    GIN_FORWARDS,
    GIN_NUM_DIRECTIONS
};

float cus_lerp(float a, float b, float t)
{
    return a + t * (b - a);
}

class BASSGinsuStream
{
private:
    // Ginsu stuff
    GinsuDataLayout GinsuHead;
    bool bDecelGin;
    uint32_t* freqSamples;
    uint32_t* cycleSamples;
    float* cycleFreqs;
    bool bLoaded;

    // BASS stuff
    HSAMPLE* hsv;
    DWORD* chans;
    HSTREAM gStream;
    
    // Playback stuff
    DWORD chan;
    DWORD channext;
    int oldcyc;
    int prevcyc;
    GIN_DIRECTION GinDir;
    float freqCurrent;
    float freqDelta;

    float SampleToFrequency(uint32_t sample)
    {
        const int segCount = GinsuHead.segCount;
        const int minFreq = GinsuHead.minFrequency;
        const int maxFreq = GinsuHead.maxFrequency;

        // Find the closest index in freqSamples
        float sampleValue = static_cast<float>(sample);
        float closestSampleValue = 0.0;
        int closestIndex = 0;
        float minDifference = std::abs(freqSamples[0] - sampleValue);

        for (int i = 1; i < segCount; ++i)
        {
            float difference = std::abs(freqSamples[i] - sampleValue);
            if (difference < minDifference)
            {
                minDifference = difference;
                closestSampleValue = freqSamples[i];
                closestIndex = i;
            }
        }

        // Calculate the frequency based on the closest interpolated value
        float startP = static_cast<float>(closestIndex) + (sampleValue - closestSampleValue) / (freqSamples[closestIndex + 1] - freqSamples[closestIndex]);
        float freq = minFreq + (maxFreq - minFreq) * (startP / segCount);

        return freq;
    }

    uint32_t FrequencyToSample(float freq)
    {
        float percentage = (freq - GinsuHead.minFrequency) / (GinsuHead.maxFrequency - GinsuHead.minFrequency);

        if (percentage < 0.0)
            percentage = 0.0;
        if (percentage > 1.0)
            percentage = 1.0;

        // Calculate the number of points in the array
        size_t numPoints = GinsuHead.segCount;

        // Calculate the index of the lower and upper points for interpolation
        float index = percentage * (numPoints);
        size_t lowerIndex = static_cast<size_t>(index);
        size_t upperIndex = std::min(lowerIndex + 1, numPoints);

        // Calculate the interpolation factor
        float interpolationFactor = index - lowerIndex;
        float interpolatedValue = cus_lerp(freqSamples[lowerIndex], freqSamples[upperIndex], interpolationFactor);

        return interpolatedValue;
    }

    uint32_t SampleToCycle(uint32_t sample, int oldCycle, bool bIsReverse)
    {
        if (bDecelGin)
            bIsReverse = !bIsReverse;

        // check for single-step first before going into a loop pass
        if (oldCycle > 0)
        {
            if (bIsReverse)
            {
                if (sample >= cycleSamples[oldCycle])
                    return oldCycle;

                if (sample >= cycleSamples[oldCycle - 1])
                {
                    return oldCycle - 1;
                }
            }
            else if (oldCycle <= (GinsuHead.cycleCount - 1))
            {
                //if (sample >= cycleSamples[oldCycle])
                //    return oldCycle;

                if (sample >= cycleSamples[oldCycle + 1])
                    return oldCycle + 1;
            }
        }

        if (sample >= cycleSamples[GinsuHead.cycleCount - 1])
            return GinsuHead.cycleCount - 1;

        for (int i = 0; i < GinsuHead.cycleCount; i++)
        {
            if (sample >= cycleSamples[GinsuHead.cycleCount - i - 1])
                return GinsuHead.cycleCount - i - 1;
        }

        return 0;
    }

    void SetFreqOnChannels(DWORD channel1, DWORD channel2, float inFreq, int cyc, bool bIsFalling)
    {
        GIN_DIRECTION loc_dir = GinDir;
        if (cyc == 0)
            loc_dir = GIN_DIRECTION::GIN_FORWARDS;
        else if (cyc >= (GinsuHead.cycleCount - 2))
            loc_dir = GIN_DIRECTION::GIN_BACKWARDS;

        float chanfreq = cycleFreqs[cyc];
        float channextfreq = cycleFreqs[cyc + (1 * loc_dir)];

        if (bIsFalling)
        {
            float t = chanfreq;
            chanfreq = channextfreq;
            channextfreq = t;
        }

        float freqscalar = inFreq / chanfreq;

        float len1 = cycleSamples[cyc + (1 * loc_dir)] - cycleSamples[cyc];
        float len22 = cycleSamples[cyc + (2 * loc_dir)] - cycleSamples[cyc + (1 * loc_dir)];

        if (bIsFalling)
        {
            double t = len1;
            len1 = len22;
            len22 = t;
        }

        float diff_fac = len22 / len1;

        float rate1 = (float)GinsuHead.sampleRate * freqscalar;
        float rate2 = (float)GinsuHead.sampleRate * diff_fac * freqscalar;
        BASS_ChannelSetAttribute(channel1, BASS_ATTRIB_FREQ, rate1);
        BASS_ChannelSetAttribute(channel2, BASS_ATTRIB_FREQ, rate2);
    }

    void SetVolOnChannels(DWORD channel1, DWORD channel2, uint32_t samp, int cyc, bool bIsFalling)
    {
        GIN_DIRECTION loc_dir = GinDir;
        if (cyc == 0)
            loc_dir = GIN_DIRECTION::GIN_FORWARDS;
        else if (cyc >= (GinsuHead.cycleCount - 2))
            loc_dir = GIN_DIRECTION::GIN_BACKWARDS;

        uint32_t ls = cycleSamples[cyc];
        uint32_t le = cycleSamples[cyc + (1 * loc_dir)];

        if (bIsFalling)
        {
            ls = cycleSamples[cyc + (1 * loc_dir)];
            le = cycleSamples[cyc];
        }

        double chan1vol = 1.0f - (((double)samp - (double)ls) / ((double)le - (double)ls));
        double chan2vol = 1.0f - chan1vol;

        BASS_ChannelSetAttribute(channel1, BASS_ATTRIB_VOL, chan1vol);
        BASS_ChannelSetAttribute(channel2, BASS_ATTRIB_VOL, chan2vol);
    }

    void SetPlaybackFrequency(float inFreq)
    {
        uint32_t samp = FrequencyToSample(inFreq);
        float cycle = SampleToCycle(samp, oldcyc, freqDelta < 0);
        int cycint = cycle;

        GIN_DIRECTION loc_dir = GinDir;
        if (cycint == 0)
            loc_dir = GIN_DIRECTION::GIN_FORWARDS;
        else if (cycint > (GinsuHead.cycleCount - 2))
            loc_dir = GIN_DIRECTION::GIN_BACKWARDS;

        // switch channels on new cycle
        if (oldcyc != cycint)
        {
            if (oldcyc < 0)
            {
                oldcyc = 0;
            }

            BASS_ChannelSetAttribute(chan, BASS_ATTRIB_VOL, 0.0f);
            BASS_ChannelSetAttribute(channext, BASS_ATTRIB_VOL, 0.0f);

            if (cycint < oldcyc)
            {
                chan = chans[cycint + (1 * loc_dir)];
                channext = chans[cycint];
            }
            else
            {
                chan = chans[cycint];
                channext = chans[cycint + (1 * loc_dir)];
            }

            BASS_ChannelSetAttribute(chan, BASS_ATTRIB_VOL, 1.0f);
            BASS_ChannelSetAttribute(channext, BASS_ATTRIB_VOL, 0.0f);

            QWORD ls = cycleSamples[cycint] * 2;
            QWORD le = cycleSamples[cycint + (1 * loc_dir)] * 2;

            QWORD ls2 = cycleSamples[cycint + (1 * loc_dir)] * 2;
            QWORD le2 = cycleSamples[cycint + (2 * loc_dir)] * 2;

            if (cycint < oldcyc)
            {
                QWORD t = ls;
                ls = ls2;
                ls2 = t;

                t = le;
                le = le2;
                le2 = t;
            }

            QWORD len = (le - ls);
            QWORD len2 = (le2 - ls2);

            QWORD basspos = BASS_ChannelGetPosition(chan, BASS_POS_BYTE);

            float dist = (float)basspos / (float)len;
            float newpos = cus_lerp(0, len2, dist);
            BASS_ChannelSetPosition(channext, (QWORD)(newpos), BASS_POS_BYTE);

            prevcyc = oldcyc;
            oldcyc = cycint;
        }

        SetFreqOnChannels(chan, channext, inFreq, cycint, cycint < prevcyc);
        SetVolOnChannels(chan, channext, samp, cycint, cycint < prevcyc);


        if (freqCurrent != inFreq)
        {
            freqDelta = inFreq - freqCurrent;
            freqCurrent = inFreq;
        }
    }

    // parsing/loading code

    bool ParseGinsu(std::filesystem::path filename)
    {
        std::ifstream ifile;
        ifile.open(filename, std::ios::binary);
        if (!ifile.is_open())
            return false;

        ifile.read((char*)&GinsuHead, sizeof(GinsuDataLayout));

        size_t sz = (GinsuHead.segCount + 1) * sizeof(uint32_t);
        freqSamples = (uint32_t*)malloc(sz);
        if (freqSamples == nullptr)
            return false;
        memset(freqSamples, 0, sz);

        for (int i = 0; i < GinsuHead.segCount + 1; i++)
        {
            uint32_t samp = 0;
            ifile.read((char*)&samp, sizeof(uint32_t));
            freqSamples[i] = samp;
        }

        // detect if it's a decel gin
        int decgintester = freqSamples[1] - freqSamples[0];
        bDecelGin = decgintester < 0;

        sz = (GinsuHead.cycleCount + 1) * sizeof(uint32_t);
        cycleSamples = (uint32_t*)malloc(sz);
        if (cycleSamples == nullptr)
            return false;
        memset(cycleSamples, 0, sz);

        for (int i = 0; i < GinsuHead.cycleCount; i++)
        {
            uint32_t cyc = 0;
            ifile.read((char*)&cyc, sizeof(uint32_t));
            cycleSamples[i] = cyc;
        }
        cycleSamples[GinsuHead.cycleCount] = UINT32_MAX;

        ifile.close();


        sz = (GinsuHead.cycleCount + 1) * sizeof(float);
        cycleFreqs = (float*)malloc(sz);
        if (cycleFreqs == nullptr)
            return false;
        memset(cycleFreqs, 0, sz);

        for (int i = 0; i < GinsuHead.cycleCount; i++)
        {
            float f = SampleToFrequency(cycleSamples[i]);
            cycleFreqs[i] = f;
        }
        cycleFreqs[GinsuHead.cycleCount] = cycleFreqs[GinsuHead.cycleCount - 1];
        return true;
    }

    bool decodeAudioStream(const char* inputFile, int outputFrequency, int16_t** outputBuffer, size_t* bufferSize) {
        VGMSTREAM* vgmstream = init_vgmstream(inputFile);
        if (!vgmstream) {
            return false; // Failed to initialize vgmstream
        }

        size_t totalSamples = (size_t)(vgmstream->num_samples);
        size_t totalChannels = (size_t)(vgmstream->channels);

        *bufferSize = totalSamples * totalChannels * sizeof(int16_t);
        *outputBuffer = (int16_t*)malloc(*bufferSize);

        int samples_to_do = (int)totalSamples;
        int16_t* out_buffer = *outputBuffer;

        // Decode the audio stream
        while (samples_to_do > 0) {
            int samples_done = render_vgmstream(out_buffer, samples_to_do, vgmstream);
            if (samples_done <= 0) {
                break; // Decoding finished or error occurred
            }

            out_buffer += samples_done * totalChannels;
            samples_to_do -= samples_done;
        }

        close_vgmstream(vgmstream);

        return true;
    }

    bool CreateChannels(std::filesystem::path filename)
    {
        int16_t* outputBuffer = nullptr;
        size_t bufferSize = 0;

        if (!decodeAudioStream(filename.string().c_str(), GinsuHead.sampleRate, &outputBuffer, &bufferSize))
            return false;

        hsv = (HSAMPLE*)malloc(GinsuHead.cycleCount * sizeof(HSAMPLE));
        if (hsv == NULL)
            return false;

        chans = (DWORD*)malloc(GinsuHead.cycleCount * sizeof(DWORD));
        if (chans == NULL)
            return false;

        for (int i = 0; i < GinsuHead.cycleCount; i++)
        {
            int16_t* cycSample = &outputBuffer[cycleSamples[i]];
            DWORD len = 0;
            if (i >= (GinsuHead.cycleCount - 1))
            {
                uintptr_t ep = (uintptr_t)(outputBuffer)+bufferSize;
                len = ep - (uintptr_t)(&outputBuffer[cycleSamples[i]]);
            }
            else
                len = (uintptr_t)(&outputBuffer[cycleSamples[i + 1]]) - (uintptr_t)(&outputBuffer[cycleSamples[i]]);

            HSAMPLE ns = BASS_SampleCreate(len, GinsuHead.sampleRate, 1, 1, BASS_SAMPLE_MONO | BASS_SAMPLE_LOOP);
            if (ns == 0)
                return false;

            if (BASS_SampleSetData(ns, cycSample) == FALSE)
                return false;

            hsv[i] = ns;

            DWORD ch = BASS_SampleGetChannel(ns, BASS_SAMCHAN_STREAM | BASS_STREAM_DECODE);
            if (ch == NULL)
                return false;

            BASS_ChannelSetAttribute(ch, BASS_ATTRIB_BUFFER, 0);
            BASS_ChannelSetAttribute(ch, BASS_ATTRIB_VOL, 0.0f);

            chans[i] = ch;
        }

        chan = chans[0];
        channext = chans[1];

        free(outputBuffer);
        return true;
    }

public:
    bool Load(std::filesystem::path ginPath)
    {
        // TODO: error handling
        if (!ParseGinsu(ginPath))
            return false;
        if (!CreateChannels(ginPath))
            return false;

        gStream = BASS_Mixer_StreamCreate(GinsuHead.sampleRate, 1, BASS_MIXER_NONSTOP | BASS_MIXER_NOSPEAKER | BASS_STREAM_DECODE);
        if (gStream == NULL)
            return false;
        
        for (int i = 0; i < GinsuHead.cycleCount; i++)
        {
            if (BASS_Mixer_StreamAddChannel(gStream, chans[i], BASS_MIXER_CHAN_NORAMPIN) == FALSE)
                return false;
        }

        BASS_ChannelSetAttribute(gStream, BASS_ATTRIB_BUFFER, 0);
        bLoaded = true;

        return true;
    }

    BOOL Play(BOOL restart)
    {
        if (!bLoaded)
            return FALSE;

        return BASS_ChannelPlay(gStream, restart);
    }

    BOOL Stop()
    {
        if (!bLoaded)
            return FALSE;

        return BASS_ChannelStop(gStream);
    }

    BOOL Pause()
    {
        if (!bLoaded)
            return FALSE;

        return BASS_ChannelPause(gStream);
    }

    BOOL SetVolume(float inVol)
    {
        if (!bLoaded)
            return FALSE;

        return BASS_ChannelSetAttribute(gStream, BASS_ATTRIB_VOL, inVol);
    }

    float GetVolume()
    {
        if (!bLoaded)
            return 0.0f;

        float ret = 0.0f;

        BOOL res = BASS_ChannelGetAttribute(gStream, BASS_ATTRIB_VOL, &ret);

        return ret;
    }

    HSTREAM GetStreamHandle()
    {
        if (!bLoaded)
            return 0;

        return gStream;
    }

    DWORD GetData(int16_t* outBuffer, DWORD length)
    {
        if (!bLoaded)
            return 0;

        return BASS_ChannelGetData(gStream, outBuffer, length);
    }

    void SetFrequency(float inFreq)
    {
        if (!bLoaded)
            return;

        SetPlaybackFrequency(inFreq);
    }

    void SetFrequencyFromPercentage(float inFreqPct)
    {
        if (!bLoaded)
            return;

        float f = cus_lerp(GinsuHead.minFrequency, GinsuHead.maxFrequency, inFreqPct);

        SetPlaybackFrequency(f);
    }

    float GetFrequency()
    {
        if (!bLoaded)
            return 0.0f;

        return freqCurrent;
    }

    float GetFrequencyFromPercentage(float inFreqPct)
    {
        if (!bLoaded)
            return 0.0f;

        float f = cus_lerp(GinsuHead.minFrequency, GinsuHead.maxFrequency, inFreqPct);

        return f;
    }

    float GetFrequencyPercentage()
    {
        if (!bLoaded)
            return 0.0f;

        float d = (freqCurrent - GinsuHead.minFrequency) / (GinsuHead.maxFrequency - GinsuHead.minFrequency);

        return d;
    }

    float GetFrequencyDelta()
    {
        if (!bLoaded)
            return 0.0f;

        return freqDelta;
    }

    bool bIsDecelGin()
    {
        if (!bLoaded)
            return false;

        return bDecelGin;
    }

    float GetMaxFrequency()
    {
        if (!bLoaded)
            return 0.0f;

        return GinsuHead.maxFrequency;
    }

    float GetMinFrequency()
    {
        if (!bLoaded)
            return 0.0f;

        return GinsuHead.minFrequency;
    }

    uint32_t GetSampleRate()
    {
        if (!bLoaded)
            return 0;

        return GinsuHead.sampleRate;
    }

    uint32_t GetSampleCount()
    {
        if (!bLoaded)
            return 0;

        return GinsuHead.sampleCount;
    }

    uint32_t GetCycleCount()
    {
        if (!bLoaded)
            return 0;

        return GinsuHead.cycleCount;
    }

    uint32_t GetSegCount()
    {
        if (!bLoaded)
            return 0;

        return GinsuHead.segCount;
    }

    BASSGinsuStream()
    {
        GinsuHead = { 0 };
        bDecelGin = false;
        freqSamples = nullptr;
        cycleSamples = nullptr;
        cycleFreqs = nullptr;
        hsv = nullptr;
        chans = nullptr;
        gStream = 0;
        oldcyc = -1;
        prevcyc = 0;
        chan = 0;
        channext = 0;
        GinDir = GIN_FORWARDS;
        freqDelta = 0.0f;
        freqCurrent = -1.0f;
        bLoaded = false;
    }
    ~BASSGinsuStream()
    {
        if (freqSamples)
            free(freqSamples);
        if (cycleSamples)
            free(cycleSamples);
        if (cycleFreqs)
            free(cycleFreqs);
        if (hsv)
            for (int i = 0; i < GinsuHead.cycleCount; i++)
            {
                BASS_SampleFree(hsv[i]);
            }
            free(hsv);
        if (chans)
        {
            free(chans);
        }
        if (gStream)
        {
            BASS_StreamFree(gStream);
        }
    }
};

// TODO: add inheritance from BASSGinsuStream and override stuff accordingly!
class BASSGinsuMultiStream
{
public:
    BASSGinsuStream* accelStream;
    BASSGinsuStream* decelStream;   


private:
    HSTREAM gStream;

    HSAMPLE hsIdle;
    HSAMPLE hsRedline;
    HSAMPLE hsReverseWhine;
    HSAMPLE hsForwardWhine;
    HSAMPLE hsIdleWhine;

    DWORD sampleRate;

    DWORD chIdle;
    DWORD chRedline;
    DWORD chReverseWhine;
    DWORD chForwardWhine;
    DWORD chIdleWhine;

    HFX fxACL;
    HFX fxDCL;
    BASS_DX8_PARAMEQ eqACL;
    BASS_DX8_PARAMEQ eqDCL;

    FrameLimiter::FPSLimitMode mFPSLimitMode;
    std::chrono::high_resolution_clock::time_point lastTime;
    double FPSLimit;

    float freqCurrent;
    float freqDelta;
    float freqDelta2;
    float freqOldDelta;

    float freqMin;
    float freqMax;

    float accelVol;
    float decelVol;
    float rateVol;
    float rateOldVol;
    float rateVolCurve;
    float rateEqCurve;
    float rateMinVol;
    float rateEqAmount;
    float rateEqOldAmount;
    float rateRPMTarget;
    float rateAccelXFadeRPMRatio;
    float rateAccelXFadeRPMRange;
    float rateDecelXFadeRPMRatio;
    float rateDecelXFadeRPMRange;
    float curRateRPMXFadeTarget;

    float rateAccelEqXFadeRPMRatio;
    float rateAccelEqXFadeRPMRange;
    float rateDecelEqXFadeRPMRatio;
    float rateDecelEqXFadeRPMRange;
    float curRateEqRPMXFadeTarget;

    float rateOfChange;

    float eqMinGainACL;
    float eqMinGainDCL;

    float eqMaxGainACL;
    float eqMaxGainDCL;

    float eqFrequencyACL;
    float eqFrequencyDCL;

    float eqBandwidthACL;
    float eqBandwidthDCL;

    float accelGlobalVol;
    float decelGlobalVol;

    float AccelXFadeRPMRatio;
    float AccelXFadeRPMRange;
    float DecelXFadeRPMRatio;
    float DecelXFadeRPMRange;
    float XFadeDelta;
    float curXFadeRange;

    float redlineStart;
    float redlineVol;
    float redlineGlobalVol;
    float redlinePitch;
    float redlineSampleRate;

    float idleEnd;
    float idleVol;
    float idleGlobalVol;
    float idlePitch;
    float idleSampleRate;

    float speedCurrent;
    float speedDelta;
    float speedOldDelta;

    float reverseWhineFadeRange;
    float reverseWhineMinVol;
    float reverseWhineVol;
    float reverseWhineGlobalVol;
    float reverseWhineSampleRate;
    bool bReverseWhineEnable;

    float forwardWhineFadeRange;
    float forwardWhineMinVol;
    float forwardWhineVol;
    float forwardWhineGlobalVol;
    float forwardWhineSampleRate;
    bool bForwardWhineEnable;

    float idleWhineFadeRange;
    float idleWhineMinVol;
    float idleWhineVol;
    float idleWhineGlobalVol;
    float idleWhineSampleRate;
    bool bIdleWhineEnable;

    bool bCurrentlyShifting;

    bool bAccelDirection;
    bool bOldDirection;

    bool bLoaded;

    BASSGinsuStream* CreateStream(std::filesystem::path ginPath)
    {
        BASSGinsuStream* stream = new BASSGinsuStream();
        if (stream == nullptr)
            return nullptr;

        if (!stream->Load(ginPath))
        {
            delete stream;
            return nullptr;
        }

        return stream;
    }

    float calculateNewSampleRate(float originalSampleRate, float semitonesShift) 
    {
        // Define the constant for the twelfth root of 2 (used to calculate semitone frequency ratio).
        //const float twelfthRootOf2 = powf(2.0, 1.0 / 12.0);
        constexpr float twelfthRootOf2 = 1.0594630943592952645618252949463f;

        // Calculate the frequency ratio for the specified number of semitones.
        float frequencyRatio = powf(twelfthRootOf2, semitonesShift);

        // Calculate the new sample rate.
        float newSampleRate = originalSampleRate * frequencyRatio;

        return newSampleRate;
    }

    float calculateRevWhineSemitones(float speedMPS)
    {
        constexpr float numerator = 16.0f / 3.6f + 1.0f;  // 16.0f * 0.22274f + 1.0f;
        constexpr float denominator = 32.0f / 3.6f + 1.0f; //32.0f * 0.22274f + 1.0f;
        constexpr float nd = numerator / denominator;
        //float result = -(3.0f / ::logf(nd)) * ::logf(0.22274f * speedKmh + 1.0f) - 11.0f;
        //float result = -(3.0f / ::logf(nd)) * ::logf(speedKmh / 3.6f + 1.0f) - 11.0f;
        float result = -(3.0f / ::logf(nd)) * ::logf(speedMPS + 1.0f) - 11.0f;
        return result;
    }

    float calculateIdleWhineSemitones(float speedMPS)
    {
        constexpr float b = 0.1402f;
        constexpr float numerator = 16.0f * b + 1.0f;
        constexpr float denominator = 32.0f * b + 1.0f;
        constexpr float nd = numerator / denominator;
        float result = -(3.0f / ::logf(nd)) * ::logf(b * (speedMPS * 3.6f) + 1.0f) - 22.0f;
        return result;
    }

    // TODO: this is a mathematical approximation. There should be a better way to calculate this (maybe based on actual gear ratios?)
    float calculateFwdWhineSemitones(float speedMPS)
    {
        constexpr float b = 0.164001f;
        //float numerator = -(3.0f / logf((16.0f * b + 1.0f) / (34.0f * b + 1.0f)));
        constexpr float numerator = 5.0348324403659308499533463488627f;

        float innerLog = logf(b * (speedMPS * 3.6f) + 1.0f);
        float result = numerator * innerLog - 12.0f;

        return result;
    }

    float linearToLogarithmic(float x, float b)
    {
        // Ensure that x is within the [0, 1.0] range
        if (x < 0.0)
        {
            x = 0.0;
        }
        else if (x > 1.0)
        {
            x = 1.0;
        }

        // Perform the logarithmic transformation
        // You can adjust the base to control the curve of the logarithmic scale
        //float base = 10.0; // Change this to your desired base (e.g., 2.0 for a different curve)
        return log10f(1.0 + (b - 1.0) * x) / log10f(b);
    }

    // TODO: make this a bit better - the bigger the frequency delta - the louder the sounds should play
    void SetPlaybackFrequency(float inFreq)
    {
        if (mFPSLimitMode == FrameLimiter::FPSLimitMode::FPS_REALTIME)
            while (!FrameLimiter::Sync_RT());
        else if (mFPSLimitMode == FrameLimiter::FPSLimitMode::FPS_ACCURATE)
            while (!FrameLimiter::Sync_SLP());

        if (freqCurrent != inFreq)
        {
            freqDelta = inFreq - freqCurrent;
            freqDelta2 = freqDelta;
            XFadeDelta = freqDelta;
            freqCurrent = inFreq;
        }
        else
        {
            freqDelta = 0.0f;
            XFadeDelta = 5.0f;
        }

        if (freqOldDelta != freqDelta)
            freqOldDelta = freqDelta;
        else
            XFadeDelta = freqOldDelta;

        std::chrono::high_resolution_clock::time_point currentTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<float> timeElapsed = std::chrono::duration_cast<std::chrono::duration<float>>(currentTime - lastTime);

        // Calculate the rate of change of RPM
        if (timeElapsed.count() > 0) 
        {
            rateOfChange = ::fabs(freqDelta2) / timeElapsed.count();
        }
        else 
        {
            rateOfChange = 0.0;  // Handle the case where timeElapsed is very small
        }

        if (freqDelta < 0)
        {
            bAccelDirection = false;
        }
        else if (freqDelta > 0)
        {
            bAccelDirection = true;
        }

        accelStream->SetFrequency(inFreq);
        if (decelStream)
        {
            decelStream->SetFrequency(inFreq);

            if (bOldDirection != bAccelDirection)
            {
                curXFadeRange = 0;
                rateOldVol = rateVol;
                rateEqOldAmount = rateEqAmount;

                if (bAccelDirection)
                {
                    curRateRPMXFadeTarget = freqCurrent + rateAccelXFadeRPMRange;
                    curRateEqRPMXFadeTarget = freqCurrent + rateAccelEqXFadeRPMRange;
                }
                else
                {
                    curRateRPMXFadeTarget = freqCurrent - rateDecelXFadeRPMRange;
                    curRateEqRPMXFadeTarget = freqCurrent - rateDecelEqXFadeRPMRange;
                }

                bOldDirection = bAccelDirection;
            }

            curXFadeRange += ::fabs(XFadeDelta);

            if (bAccelDirection)
            {
                if (curXFadeRange > AccelXFadeRPMRange)
                    curXFadeRange = AccelXFadeRPMRange;

                float d = curXFadeRange / AccelXFadeRPMRange;
                accelVol = d;
                decelVol = 1.0f - d;
            }
            else
            {
                if (curXFadeRange > DecelXFadeRPMRange)
                    curXFadeRange = DecelXFadeRPMRange;

                float d = curXFadeRange / DecelXFadeRPMRange;
                decelVol = d;
                accelVol = 1.0f - d;

                if ((d >= 1.0f) || (inFreq < decelStream->GetMinFrequency()))
                {
                    // as we go closer to the min frequency of decel gin, we crossfade back into accel gin
                    float dclDist = (inFreq - decelStream->GetMinFrequency()) / (decelStream->GetMaxFrequency() - decelStream->GetMinFrequency());
                    float aclDist = 1.0f - dclDist;

                    decelVol = dclDist;
                    accelVol = aclDist;
                }
            }

            decelVol = std::clamp(decelVol, 0.0f, 1.0f);
            accelVol = std::clamp(accelVol, 0.0f, 1.0f);
        }


        // TODO: make this a teeny bit better; add some accumulator that will fill up from redline start and add up quicker closer it is to the max RPM
        // as the accumulator gets closer to max, louder the redline sound is
        if (chRedline)
        {
            float d = (inFreq - redlineStart) / (freqMax - redlineStart);
            redlineVol = d;
            redlineVol = std::clamp(redlineVol, 0.0f, 1.0f);
            redlineVol *= redlineGlobalVol;

            float sensitivityFactor = 1.0f - ::powf(redlineVol, 2);

            decelVol *= sensitivityFactor;
            accelVol *= sensitivityFactor;

            decelVol = std::clamp(decelVol, 0.0f, 1.0f);
            accelVol = std::clamp(accelVol, 0.0f, 1.0f);

            BASS_ChannelSetAttribute(chRedline, BASS_ATTRIB_VOL, redlineVol);
        }

        if (chIdle)
        {
            float dist = freqMax - freqMin;
            float p = (inFreq - freqMin);
            float d = p / idleEnd;

            idleVol = 1.0f - d;
            idleVol = std::clamp(idleVol, 0.0f, 1.0f);
            idleVol *= idleGlobalVol;

            decelVol *= 1.0f - idleVol;
            accelVol *= 1.0f - idleVol;

            decelVol = std::clamp(decelVol, 0.0f, 1.0f);
            accelVol = std::clamp(accelVol, 0.0f, 1.0f);

            BASS_ChannelSetAttribute(chIdle, BASS_ATTRIB_VOL, idleVol);
        }

        float ratePct = std::clamp(rateOfChange / rateRPMTarget, 0.0f, 1.0f);
        float rateVolCalc = std::clamp(cus_lerp(rateMinVol, 1.0f, linearToLogarithmic(ratePct, rateVolCurve)), 0.0f, 1.0f);

        if (bAccelDirection)
        {
            float d = std::clamp(freqCurrent / curRateRPMXFadeTarget, 0.0f, 1.0f);
            rateVol = cus_lerp(rateOldVol, rateVolCalc, d);
        }
        else
        {
            float d = std::clamp(curRateRPMXFadeTarget / freqCurrent, 0.0f, 1.0f);
            rateVol = cus_lerp(rateOldVol, rateVolCalc, d);
        }

        float rateEqAmountCalc = 0;

        if (rateEqCurve <= 1.0f)
            rateEqAmountCalc = std::clamp(cus_lerp(0.0f, 1.0f, ratePct), 0.0f, 1.0f);
        else
            rateEqAmountCalc = std::clamp(cus_lerp(0.0f, 1.0f, linearToLogarithmic(ratePct, rateEqCurve)), 0.0f, 1.0f);
        
        if (bAccelDirection)
        {
            float d = std::clamp(freqCurrent / curRateEqRPMXFadeTarget, 0.0f, 1.0f);
            rateEqAmount = cus_lerp(rateEqOldAmount, rateEqAmountCalc, d);
        }
        else
        {
            float d = std::clamp(curRateEqRPMXFadeTarget / freqCurrent, 0.0f, 1.0f);
            rateEqAmount = cus_lerp(rateEqOldAmount, rateEqAmountCalc, d);
        }

        eqACL.fGain = cus_lerp(eqMinGainACL, eqMaxGainACL, rateEqAmount);

        if (decelStream)
        {
            decelVol *= rateVol;
            eqDCL.fGain = cus_lerp(eqMinGainDCL, eqMaxGainDCL, rateEqAmount);
            BASS_FXSetParameters(fxDCL, &eqDCL);
            decelStream->SetVolume(decelVol * decelGlobalVol);
        }

        accelVol *= rateVol;

        lastTime = currentTime;

        BASS_FXSetParameters(fxACL, &eqACL);
        accelStream->SetVolume(accelVol * accelGlobalVol);
    }

    void SetPlaybackSpeed(float inSpeedMPS)
    {
        if (!bForwardWhineEnable && !bReverseWhineEnable && !bIdleWhineEnable)
            return;

        if (!chForwardWhine && !chReverseWhine && !chIdleWhine)
            return;

        if (speedCurrent != inSpeedMPS)
        {
            speedDelta = inSpeedMPS - speedCurrent;
            speedCurrent = inSpeedMPS;
        }
        else
            speedDelta = 0.0f;

        if (speedOldDelta != speedDelta)
            speedOldDelta = speedDelta;

        if (chForwardWhine)
        {
            if (bForwardWhineEnable && !bCurrentlyShifting)
            {
                forwardWhineVol = std::clamp(cus_lerp(forwardWhineMinVol, 1.0f, inSpeedMPS / forwardWhineFadeRange), 0.0f, 1.0f);
                

                BASS_ChannelSetAttribute(chForwardWhine, BASS_ATTRIB_VOL, forwardWhineVol * forwardWhineGlobalVol);
                BASS_ChannelSetAttribute(chForwardWhine, BASS_ATTRIB_FREQ, calculateNewSampleRate(forwardWhineSampleRate, calculateFwdWhineSemitones(inSpeedMPS)));
            }
            else
                BASS_ChannelSetAttribute(chForwardWhine, BASS_ATTRIB_VOL, 0.0f);
        }

        if (chIdleWhine)
        {
            if (bIdleWhineEnable && !bCurrentlyShifting)
            {
                idleWhineVol = 1.0f - std::clamp(cus_lerp(idleWhineMinVol, 1.0f, inSpeedMPS / idleWhineFadeRange), 0.0f, 1.0f);


                BASS_ChannelSetAttribute(chIdleWhine, BASS_ATTRIB_VOL, idleWhineVol * idleWhineGlobalVol);
                BASS_ChannelSetAttribute(chIdleWhine, BASS_ATTRIB_FREQ, calculateNewSampleRate(idleWhineSampleRate, calculateIdleWhineSemitones(inSpeedMPS)));
            }
            else
                BASS_ChannelSetAttribute(chIdleWhine, BASS_ATTRIB_VOL, 0.0f);
        }

        if (chReverseWhine)
        {
            if (bReverseWhineEnable && !bCurrentlyShifting)
            {
                reverseWhineVol = std::clamp(cus_lerp(reverseWhineMinVol, 1.0f, inSpeedMPS / reverseWhineFadeRange), 0.0f, 1.0f);

                BASS_ChannelSetAttribute(chReverseWhine, BASS_ATTRIB_VOL, reverseWhineVol * reverseWhineGlobalVol);
                BASS_ChannelSetAttribute(chReverseWhine, BASS_ATTRIB_FREQ, calculateNewSampleRate(reverseWhineSampleRate, calculateRevWhineSemitones(inSpeedMPS)));
            }
            else
                BASS_ChannelSetAttribute(chReverseWhine, BASS_ATTRIB_VOL, 0.0f);
        }
    }

    void ReleaseEverything()
    {
        bLoaded = false;
        if (gStream)
            BASS_StreamFree(gStream);
        if (chRedline)
            BASS_ChannelFree(chRedline);
        if (hsRedline)
            BASS_SampleFree(hsRedline);
        if (chIdle)
            BASS_ChannelFree(chIdle);
        if (hsIdle)
            BASS_SampleFree(hsIdle);
        if (chReverseWhine)
            BASS_ChannelFree(chReverseWhine);
        if (hsReverseWhine)
            BASS_SampleFree(hsReverseWhine);
        if (chForwardWhine)
            BASS_ChannelFree(chForwardWhine);
        if (hsForwardWhine)
            BASS_SampleFree(hsForwardWhine);
        if (chIdleWhine)
            BASS_ChannelFree(chIdleWhine);
        if (hsIdleWhine)
            BASS_SampleFree(hsIdleWhine);
        if (accelStream)
            delete accelStream;
        if (decelStream)
            delete decelStream;
    }

    void SetDefaultRedline()
    {
        // set the redline by default to be 5 percent of the max RPM
        redlineStart = freqMax - (freqMax * 0.05f);
    }

    void SetDefaultIdle()
    {
        // set the idle end by default to be 7.7 percent of the total distance
        idleEnd = (freqMax - freqMin) * 0.077f;
    }

public:
    // Loads the sounds and creates the streams.
    // Should only be called from BASSGinsuPlayer
    bool Load(DWORD freq, std::filesystem::path accelGinPath, std::filesystem::path decelGinPath, std::filesystem::path redlinePath, std::filesystem::path idlePath, std::filesystem::path reverseWhinePath, std::filesystem::path forwardWhinePath, std::filesystem::path idleWhinePath)
    {
        sampleRate = freq;

        accelStream = CreateStream(accelGinPath);
        if (accelStream == nullptr)
            return false;

        if (!decelGinPath.empty())
        {
            decelStream = CreateStream(decelGinPath);
            if (decelStream == nullptr)
            {
                if (accelStream)
                    delete accelStream;
                return false;
            }
        }

        if (!redlinePath.empty())
        {
            hsRedline = BASS_SampleLoad(FALSE, redlinePath.u16string().c_str(), 0, 0, 1, BASS_SAMPLE_LOOP | BASS_SAMPLE_MONO | BASS_UNICODE);
            if (hsRedline == 0)
            {
                ReleaseEverything();
                return false;
            }
            chRedline = BASS_SampleGetChannel(hsRedline, BASS_SAMCHAN_STREAM | BASS_STREAM_DECODE);
            if (chRedline == NULL)
            {
                ReleaseEverything();
                return false;
            }
            BASS_ChannelSetAttribute(chRedline, BASS_ATTRIB_VOL, redlineVol * redlineGlobalVol);
        }

        if (!idlePath.empty())
        {
            hsIdle = BASS_SampleLoad(FALSE, idlePath.u16string().c_str(), 0, 0, 1, BASS_SAMPLE_LOOP | BASS_SAMPLE_MONO | BASS_UNICODE);
            if (hsIdle == 0)
            {
                ReleaseEverything();
                return false;
            }
            chIdle = BASS_SampleGetChannel(hsIdle, BASS_SAMCHAN_STREAM | BASS_STREAM_DECODE);
            if (chIdle == NULL)
            {
                ReleaseEverything();
                return false;
            }
            BASS_ChannelSetAttribute(chIdle, BASS_ATTRIB_VOL, idleVol * idleGlobalVol);
        }

        if (!reverseWhinePath.empty())
        {
            hsReverseWhine = BASS_SampleLoad(FALSE, reverseWhinePath.u16string().c_str(), 0, 0, 1, BASS_SAMPLE_LOOP | BASS_SAMPLE_MONO | BASS_UNICODE);
            if (hsReverseWhine == 0)
            {
                ReleaseEverything();
                return false;
            }
            chReverseWhine = BASS_SampleGetChannel(hsReverseWhine, BASS_SAMCHAN_STREAM | BASS_STREAM_DECODE);
            if (chReverseWhine == NULL)
            {
                ReleaseEverything();
                return false;
            }
            BASS_ChannelSetAttribute(chReverseWhine, BASS_ATTRIB_VOL, 0.0f);
            BASS_ChannelPause(chReverseWhine);
        }

        if (!forwardWhinePath.empty())
        {
            hsForwardWhine = BASS_SampleLoad(FALSE, forwardWhinePath.u16string().c_str(), 0, 0, 1, BASS_SAMPLE_LOOP | BASS_SAMPLE_MONO | BASS_UNICODE);
            if (hsForwardWhine == 0)
            {
                ReleaseEverything();
                return false;
            }
            chForwardWhine = BASS_SampleGetChannel(hsForwardWhine, BASS_SAMCHAN_STREAM | BASS_STREAM_DECODE);
            if (chForwardWhine == NULL)
            {
                ReleaseEverything();
                return false;
            }
            BASS_ChannelSetAttribute(chForwardWhine, BASS_ATTRIB_VOL, 0.0f);
            BASS_ChannelPause(chForwardWhine);
        }

        if (!idleWhinePath.empty())
        {
            hsIdleWhine = BASS_SampleLoad(FALSE, idleWhinePath.u16string().c_str(), 0, 0, 1, BASS_SAMPLE_LOOP | BASS_SAMPLE_MONO | BASS_UNICODE);
            if (hsIdleWhine == 0)
            {
                ReleaseEverything();
                return false;
            }
            chIdleWhine = BASS_SampleGetChannel(hsIdleWhine, BASS_SAMCHAN_STREAM | BASS_STREAM_DECODE);
            if (chIdleWhine == NULL)
            {
                ReleaseEverything();
                return false;
            }
            BASS_ChannelSetAttribute(chIdleWhine, BASS_ATTRIB_VOL, 0.0f);
            BASS_ChannelPause(chIdleWhine);
        }

        // we'll use accelStream's sample rate as a basis -- TODO: push it further if one of the other samples need it
        gStream = BASS_Mixer_StreamCreate(sampleRate, 1, BASS_MIXER_NONSTOP | BASS_MIXER_NOSPEAKER | BASS_STREAM_DECODE);
        if (gStream == NULL)
        {
            ReleaseEverything();
            return false;
        }

        // plug the streams together
        if (BASS_Mixer_StreamAddChannel(gStream, accelStream->GetStreamHandle(), BASS_MIXER_CHAN_NORAMPIN) == FALSE)
        {
            ReleaseEverything();
            return false;
        }

        if (decelStream)
        {
            if (BASS_Mixer_StreamAddChannel(gStream, decelStream->GetStreamHandle(), BASS_MIXER_CHAN_NORAMPIN) == FALSE)
            {
                ReleaseEverything();
                return false;
            }
        }

        if (chRedline)
        {
            if (BASS_Mixer_StreamAddChannel(gStream, chRedline, BASS_MIXER_CHAN_NORAMPIN) == FALSE)
            {
                ReleaseEverything();
                return false;
            }
        }

        if (chIdle)
        {
            if (BASS_Mixer_StreamAddChannel(gStream, chIdle, BASS_MIXER_CHAN_NORAMPIN) == FALSE)
            {
                ReleaseEverything();
                return false;
            }
        }

        if (chReverseWhine)
        {
            if (BASS_Mixer_StreamAddChannel(gStream, chReverseWhine, BASS_MIXER_CHAN_NORAMPIN) == FALSE)
            {
                ReleaseEverything();
                return false;
            }
        }

        if (chForwardWhine)
        {
            if (BASS_Mixer_StreamAddChannel(gStream, chForwardWhine, BASS_MIXER_CHAN_NORAMPIN) == FALSE)
            {
                ReleaseEverything();
                return false;
            }
        }

        if (chIdleWhine)
        {
            if (BASS_Mixer_StreamAddChannel(gStream, chIdleWhine, BASS_MIXER_CHAN_NORAMPIN) == FALSE)
            {
                ReleaseEverything();
                return false;
            }
        }

        // accelStream takes precedence, so decelStream is silent first
        accelStream->SetVolume(accelVol);
        if (decelStream)
            decelStream->SetVolume(decelVol);

        BASS_ChannelSetAttribute(gStream, BASS_ATTRIB_BUFFER, 0);

        freqMin = accelStream->GetMinFrequency();
        freqMax = accelStream->GetMaxFrequency();

        if (decelStream)
        {
            SetAccelXFadeRatio(AccelXFadeRPMRatio);
            SetDecelXFadeRatio(DecelXFadeRPMRatio);
            SetRateAccelXFadeRatio(rateAccelXFadeRPMRatio);
            SetRateDecelXFadeRatio(rateDecelXFadeRPMRatio);
            SetRateAccelEqXFadeRatio(rateAccelEqXFadeRPMRatio);
            SetRateDecelEqXFadeRatio(rateDecelEqXFadeRPMRatio);
        }

        if (chRedline)
        {
            BASS_ChannelGetAttribute(chRedline, BASS_ATTRIB_FREQ, &redlineSampleRate);
            SetDefaultRedline();
        }

        if (chIdle)
        {
            BASS_ChannelGetAttribute(chIdle, BASS_ATTRIB_FREQ, &idleSampleRate);
            SetDefaultIdle();
        }

        if (chReverseWhine)
        {
            BASS_ChannelGetAttribute(chReverseWhine, BASS_ATTRIB_FREQ, &reverseWhineSampleRate);
        }

        if (chForwardWhine)
        {
            BASS_ChannelGetAttribute(chForwardWhine, BASS_ATTRIB_FREQ, &forwardWhineSampleRate);
        }

        if (chIdleWhine)
        {
            BASS_ChannelGetAttribute(chIdleWhine, BASS_ATTRIB_FREQ, &idleWhineSampleRate);
        }

        // set up BASS_FX
        fxACL = BASS_ChannelSetFX(accelStream->GetStreamHandle(), BASS_FX_DX8_PARAMEQ, 0);
        if (fxACL)
        {
            BASS_FXSetParameters(fxACL, &eqACL);
        }
        
        if (decelStream)
        {
            fxDCL = BASS_ChannelSetFX(decelStream->GetStreamHandle(), BASS_FX_DX8_PARAMEQ, 0);
            if (fxDCL)
            {
                BASS_FXSetParameters(fxDCL, &eqDCL);
            }
        }

        // set up FrameLimiter
        FrameLimiter::Init(mFPSLimitMode, FPSLimit);

        bLoaded = true;
        return true;
    }

    BOOL Play(BOOL restart)
    {
        if (!bLoaded)
            return FALSE;

        return BASS_ChannelPlay(gStream, restart);
    }

    BOOL Stop()
    {
        if (!bLoaded)
            return FALSE;

        return BASS_ChannelStop(gStream);
    }

    BOOL Pause()
    {
        if (!bLoaded)
            return FALSE;

        return BASS_ChannelPause(gStream);
    }

    BOOL SetVolume(float inVol)
    {
        if (!bLoaded)
            return FALSE;

        return BASS_ChannelSetAttribute(gStream, BASS_ATTRIB_VOL, inVol);
    }

    float GetVolume()
    {
        if (!bLoaded)
            return 0.0f;

        float ret = 0.0f;

        BOOL res = BASS_ChannelGetAttribute(gStream, BASS_ATTRIB_VOL, &ret);

        return ret;
    }

    // Set the volume of the accel gin stream
    void SetAccelVolume(float inVol)
    {
        if (!bLoaded)
            return;

        accelGlobalVol = std::clamp(inVol, 0.0f, 1.0f);
    }

    // Set the volume of the decel gin stream
    void SetDecelVolume(float inVol)
    {
        if (!bLoaded)
            return;
        if (!decelStream)
            return;

        decelGlobalVol = std::clamp(inVol, 0.0f, 1.0f);
    }

    // Set the volume of the redline sound
    void SetRedlineVolume(float inVol)
    {
        if (!bLoaded)
            return;
        if (!chRedline)
            return;

        redlineGlobalVol = std::clamp(inVol, 0.0f, 1.0f);
    }

    // Set the volume of the idle sound
    void SetIdleVolume(float inVol)
    {
        if (!bLoaded)
            return;
        if (!chIdle)
            return;

        idleGlobalVol = std::clamp(inVol, 0.0f, 1.0f);
    }

    void SetReverseWhineVolume(float inVol)
    {
        if (!bLoaded)
            return;
        if (!chIdle)
            return;

        reverseWhineGlobalVol = std::clamp(inVol, 0.0f, 1.0f);
    }

    void SetForwardWhineVolume(float inVol)
    {
        if (!bLoaded)
            return;
        if (!chIdle)
            return;

        forwardWhineGlobalVol = std::clamp(inVol, 0.0f, 1.0f);
    }

    void SetIdleWhineVolume(float inVol)
    {
        if (!bLoaded)
            return;
        if (!chIdle)
            return;

        idleWhineGlobalVol = std::clamp(inVol, 0.0f, 1.0f);
    }

    HSTREAM GetStreamHandle()
    {
        if (!bLoaded)
            return 0;

        return gStream;
    }

    // Gets the data from the stream into the outBuffer
    DWORD GetData(int16_t* outBuffer, DWORD length)
    {
        if (!bLoaded)
            return 0;

        return BASS_ChannelGetData(gStream, outBuffer, length);
    }

    void SetFrequency(float inFreq)
    {
        if (!bLoaded)
            return;

        SetPlaybackFrequency(inFreq);
    }

    float GetFrequency()
    {
        if (!bLoaded)
            return 0.0f;

        return freqCurrent;
    }

    float GetFrequencyDelta()
    {
        if (!bLoaded)
            return 0.0f;

        return freqDelta;
    }

    // Gets the maximum RPM
    float GetMaxFrequency()
    {
        if (!bLoaded)
            return 0.0f;

        return freqMax;
    }

    // Gets the minimum RPM
    float GetMinFrequency()
    {
        if (!bLoaded)
            return 0.0f;

        return freqMin;
    }

    // Sets the maximum RPM 
    void SetMaxFrequency(float inFreq)
    {
        if (!bLoaded)
            return;

        freqMax = inFreq;

        if (decelStream)
        {
            SetAccelXFadeRatio(AccelXFadeRPMRatio);
            SetDecelXFadeRatio(DecelXFadeRPMRatio);
            SetRateAccelXFadeRatio(rateAccelXFadeRPMRatio);
            SetRateDecelXFadeRatio(rateDecelXFadeRPMRatio);
            SetRateAccelEqXFadeRatio(rateAccelEqXFadeRPMRatio);
            SetRateDecelEqXFadeRatio(rateDecelEqXFadeRPMRatio);
        }

        if (chRedline)
            SetDefaultRedline();

        if (chIdle)
            SetDefaultIdle();

        return;
    }

    // Sets the minimum RPM 
    void SetMinFrequency(float inFreq)
    {
        if (!bLoaded)
            return;

        freqMin = inFreq;

        if (decelStream)
        {
            SetAccelXFadeRatio(AccelXFadeRPMRatio);
            SetDecelXFadeRatio(DecelXFadeRPMRatio);
            SetRateAccelXFadeRatio(rateAccelXFadeRPMRatio);
            SetRateDecelXFadeRatio(rateDecelXFadeRPMRatio);
            SetRateAccelEqXFadeRatio(rateAccelEqXFadeRPMRatio);
            SetRateDecelEqXFadeRatio(rateDecelEqXFadeRPMRatio);
        }

        if (chRedline)
            SetDefaultRedline();

        if (chIdle)
            SetDefaultIdle();

        return;
    }

    void SetSpeed(float inSpeedMPS)
    {
        if (!bLoaded)
            return;

        SetPlaybackSpeed(inSpeedMPS);
    }

    float GetSpeed()
    {
        if (!bLoaded)
            return 0.0f;

        return speedCurrent;
    }

    void SetCurrentlyShifting(bool bStatus)
    {
        if (!bLoaded)
            return;

        bCurrentlyShifting = bStatus;
    }

    bool GetCurrentlyShifting()
    {
        if (!bLoaded)
            return false;

        return bCurrentlyShifting;
    }

    // Sets at which RPM the redline sound should start fading in
    void SetRedlineStart(float inFreq)
    {
        if (!bLoaded)
            return;

        redlineStart = inFreq;

        return;
    }

    // Sets at which RPM the idle sound should start fading out
    void SetIdleEnd(float inFreq)
    {
        if (!bLoaded)
            return;

        idleEnd = inFreq;

        return;
    }

    BOOL SetRedlinePitch(float inPitch)
    {
        if (!bLoaded)
            return FALSE;
        if (!chRedline)
            return FALSE;

        redlinePitch = inPitch;

        return BASS_ChannelSetAttribute(chRedline, BASS_ATTRIB_FREQ, redlineSampleRate * redlinePitch);
    }

    BOOL SetIdlePitch(float inPitch)
    {
        if (!bLoaded)
            return FALSE;
        if (!chIdle)
            return FALSE;

        idlePitch = inPitch;

        return BASS_ChannelSetAttribute(chIdle, BASS_ATTRIB_FREQ, idleSampleRate * idlePitch);
    }

    uint32_t GetSampleRate()
    {
        if (!bLoaded)
            return 0;

        float v = 0.0f;
        BASS_ChannelGetAttribute(gStream, BASS_ATTRIB_FREQ, &v);
        uint32_t sampleRate = (uint32_t)v;

        return sampleRate;
    }

    // Sets the decel crossfade range based on the RPM range divided by inRatio
    // This is the "buffer" range between accel and decel sound transitions
    void SetDecelXFadeRatio(float inRatio)
    {
        if (!bLoaded)
            return;

        if (!decelStream)
            return;

        DecelXFadeRPMRatio = inRatio;

        float RPMrange = freqMax - freqMin;
        DecelXFadeRPMRange = RPMrange / DecelXFadeRPMRatio;
    }

    // Gets the decel crossfade range ratio
    float GetDecelXFadeRatio()
    {
        if (!bLoaded)
            return 0.0f;

        if (!decelStream)
            return 0.0f;

        return DecelXFadeRPMRatio;
    }

    // Sets the decel crossfade range based on the RPM range divided by inRatio
    // This is the "buffer" range between accel and decel sound transitions
    void SetAccelXFadeRatio(float inRatio)
    {
        if (!bLoaded)
            return;

        if (!decelStream)
            return;

        AccelXFadeRPMRatio = inRatio;

        float RPMrange = freqMax - freqMin;
        AccelXFadeRPMRange = RPMrange / AccelXFadeRPMRatio;
    }

    // Gets the decel crossfade range ratio
    float GetAccelXFadeRatio()
    {
        if (!bLoaded)
            return 0.0f;

        if (!decelStream)
            return 0.0f;

        return AccelXFadeRPMRatio;
    }

    // Sets the decel crossfade range manually
    // This is the "buffer" range between accel and decel sound transitions
    void SetDecelXFadeRange(float inRange)
    {
        if (!bLoaded)
            return;

        DecelXFadeRPMRange = inRange;
    }

    // Gets the current decel crossfade range
    float GetDecelXFadeRange()
    {
        if (!bLoaded)
            return 0.0f;

        if (!decelStream)
            return 0.0f;

        return DecelXFadeRPMRange;
    }

    // Sets the accel crossfade range manually
    // This is the "buffer" range between accel and decel sound transitions
    void SetAccelXFadeRange(float inRange)
    {
        if (!bLoaded)
            return;

        AccelXFadeRPMRange = inRange;
    }

    // Gets the current accel crossfade range
    float GetAccelXFadeRange()
    {
        if (!bLoaded)
            return 0.0f;

        if (!decelStream)
            return 0.0f;

        return AccelXFadeRPMRange;
    }
    // Sets the decel crossfade range based on the RPM range divided by inRatio
    // This is the "buffer" range between accel and decel sound transitions
    void SetRateDecelXFadeRatio(float inRatio)
    {
        if (!bLoaded)
            return;

        if (!decelStream)
            return;

        rateDecelXFadeRPMRatio = inRatio;

        float RPMrange = freqMax - freqMin;
        rateDecelXFadeRPMRange = RPMrange / rateDecelXFadeRPMRatio;
    }

    // Gets the decel crossfade range ratio
    float GetRateDecelXFadeRatio()
    {
        if (!bLoaded)
            return 0.0f;

        if (!decelStream)
            return 0.0f;

        return rateDecelXFadeRPMRatio;
    }

    // Sets the decel crossfade range based on the RPM range divided by inRatio
    // This is the "buffer" range between accel and decel sound transitions
    void SetRateAccelXFadeRatio(float inRatio)
    {
        if (!bLoaded)
            return;

        if (!decelStream)
            return;

        rateAccelXFadeRPMRatio = inRatio;

        float RPMrange = freqMax - freqMin;
        rateAccelXFadeRPMRange = RPMrange / rateAccelXFadeRPMRatio;
    }

    // Gets the decel crossfade range ratio
    float GetRateAccelXFadeRatio()
    {
        if (!bLoaded)
            return 0.0f;

        if (!decelStream)
            return 0.0f;

        return rateAccelXFadeRPMRatio;
    }

    // Sets the decel crossfade range manually
    // This is the "buffer" range between accel and decel sound transitions
    void SetRateDecelXFadeRange(float inRange)
    {
        if (!bLoaded)
            return;

        rateDecelXFadeRPMRange = inRange;
    }

    float GetRateDecelXFadeRange()
    {
        if (!bLoaded)
            return 0.0f;

        if (!decelStream)
            return 0.0f;

        return rateDecelXFadeRPMRange;
    }

    void SetRateAccelXFadeRange(float inRange)
    {
        if (!bLoaded)
            return;

        rateAccelXFadeRPMRange = inRange;
    }

    float GetRateAccelXFadeRange()
    {
        if (!bLoaded)
            return 0.0f;

        if (!decelStream)
            return 0.0f;

        return rateAccelXFadeRPMRange;
    }

    void SetRateDecelEqXFadeRatio(float inRatio)
    {
        if (!bLoaded)
            return;

        if (!decelStream)
            return;

        rateDecelEqXFadeRPMRatio = inRatio;

        float RPMrange = freqMax - freqMin;
        rateDecelEqXFadeRPMRange = RPMrange / rateDecelEqXFadeRPMRatio;
    }

    // Gets the decel crossfade range ratio
    float GetRateDecelEqXFadeRatio()
    {
        if (!bLoaded)
            return 0.0f;

        if (!decelStream)
            return 0.0f;

        return rateDecelEqXFadeRPMRatio;
    }

    // Sets the decel crossfade range based on the RPM range divided by inRatio
    // This is the "buffer" range between accel and decel sound transitions
    void SetRateAccelEqXFadeRatio(float inRatio)
    {
        if (!bLoaded)
            return;

        if (!decelStream)
            return;

        rateAccelEqXFadeRPMRatio = inRatio;

        float RPMrange = freqMax - freqMin;
        rateAccelEqXFadeRPMRange = RPMrange / rateAccelEqXFadeRPMRatio;
    }

    // Gets the decel crossfade range ratio
    float GetRateAccelEqXFadeRatio()
    {
        if (!bLoaded)
            return 0.0f;

        if (!decelStream)
            return 0.0f;

        return rateAccelEqXFadeRPMRatio;
    }


    void SetRateDecelEqXFadeRange(float inRange)
    {
        if (!bLoaded)
            return;

        rateDecelEqXFadeRPMRange = inRange;
    }

    float GetRateDecelEqXFadeRange()
    {
        if (!bLoaded)
            return 0.0f;

        if (!decelStream)
            return 0.0f;

        return rateDecelEqXFadeRPMRange;
    }

    void SetRateAccelEqXFadeRange(float inRange)
    {
        if (!bLoaded)
            return;

        rateAccelEqXFadeRPMRange = inRange;
    }

    float GetRateAccelEqXFadeRange()
    {
        if (!bLoaded)
            return 0.0f;

        if (!decelStream)
            return 0.0f;

        return rateAccelEqXFadeRPMRange;
    }

    void SetReverseWhineEnable(bool bEnable)
    {
        if (!bLoaded)
            return;
        if (chReverseWhine)
        {
            bReverseWhineEnable = bEnable;
            if (bEnable)
            {
                BASS_ChannelPlay(chReverseWhine, FALSE);
            }
            else
            {
                BASS_ChannelSetAttribute(chReverseWhine, BASS_ATTRIB_VOL, 0.0f);
                BASS_ChannelPause(chReverseWhine);
            }
        }
    }

    void SetForwardWhineEnable(bool bEnable)
    {
        if (!bLoaded)
            return;
        if (chForwardWhine)
        {
            bForwardWhineEnable = bEnable;
            if (bEnable)
            {
                BASS_ChannelPlay(chForwardWhine, FALSE);
            }
            else
            {
                BASS_ChannelSetAttribute(chForwardWhine, BASS_ATTRIB_VOL, 0.0f);
                BASS_ChannelPause(chForwardWhine);
            }
        }
    }

    void SetIdleWhineEnable(bool bEnable)
    {
        if (!bLoaded)
            return;
        if (chIdleWhine)
        {
            bIdleWhineEnable = bEnable;
            if (bEnable)
            {
                BASS_ChannelPlay(chIdleWhine, FALSE);
            }
            else
            {
                BASS_ChannelSetAttribute(chIdleWhine, BASS_ATTRIB_VOL, 0.0f);
                BASS_ChannelPause(chIdleWhine);
            }
        }
    }

    bool GetReverseWhineEnable()
    {
        if (!bLoaded)
            return false;
        if (!chReverseWhine)
            return false;

        return bReverseWhineEnable;
    }

    bool GetForwardWhineEnable()
    {
        if (!bLoaded)
            return false;
        if (!chForwardWhine)
            return false;

        return bForwardWhineEnable;
    }

    bool GetIdleWhineEnable()
    {
        if (!bLoaded)
            return false;
        if (!chIdleWhine)
            return false;

        return bIdleWhineEnable;
    }

    void SetReverseWhineMinVol(float inVol)
    {
        if (!bLoaded)
            return;
        if (chReverseWhine)
            reverseWhineMinVol = inVol;
    }

    void SetForwardWhineMinVol(float inVol)
    {
        if (!bLoaded)
            return;
        if (chForwardWhine)
            forwardWhineMinVol = inVol;
    }

    void SetIdleWhineMinVol(float inVol)
    {
        if (!bLoaded)
            return;
        if (chIdleWhine)
            idleWhineMinVol = inVol;
    }

    float GetIdleWhineMinVol()
    {
        if (!bLoaded)
            return 0.0f;

        if (!chIdleWhine)
            return 0.0f;

        return idleWhineMinVol;
    }

    void SetReverseWhineFadeRange(float inSpeedMPS)
    {
        if (!bLoaded)
            return;
        if (chReverseWhine)
            reverseWhineFadeRange = inSpeedMPS;
    }

    void SetIdleWhineFadeRange(float inSpeedMPS)
    {
        if (!bLoaded)
            return;
        if (chIdleWhine)
            idleWhineFadeRange = inSpeedMPS;
    }

    float GetReverseWhineFadeRange()
    {
        if (!bLoaded)
            return 0.0f;
        if (!chReverseWhine)
            return 0.0f;
        return reverseWhineFadeRange;
    }

    float GetIdleWhineFadeRange()
    {
        if (!bLoaded)
            return 0.0f;
        if (!chIdleWhine)
            return 0.0f;
        return idleWhineFadeRange;
    }

    void SetRateMinVol(float inVol)
    {
        if (!bLoaded)
            return;
        rateMinVol = inVol;
    }

    float GetRateMinVol()
    {
        if (!bLoaded)
            return 0.0f;
        return rateMinVol;
    }

    void SetRateEqCurve(float inCurve)
    {
        if (!bLoaded)
            return;
        rateEqCurve = inCurve;
    }

    float GetRateEqCurve()
    {
        if (!bLoaded)
            return 0.0f;
        return rateEqCurve;
    }

    void SetRateVolCurve(float inCurve)
    {
        if (!bLoaded)
            return;
        rateVolCurve = inCurve;
    }

    float GetRateVolCurve()
    {
        if (!bLoaded)
            return 0.0f;
        return rateVolCurve;
    }

    void SetRateRPMTarget(float inRPM)
    {
        if (!bLoaded)
            return;
        rateRPMTarget = inRPM;
    }

    float GetRateRPMTarget()
    {
        if (!bLoaded)
            return 0.0f;
        return rateRPMTarget;
    }

    void SetFPSLimit(double inFPS)
    {
        if (!bLoaded)
            return;
        FPSLimit = inFPS;
    }

    double GetFPSLimit()
    {
        if (!bLoaded)
            return 0.0f;
        return FPSLimit;
    }

    BASSGinsuMultiStream()
    {
        mFPSLimitMode = FrameLimiter::FPSLimitMode::FPS_ACCURATE;

        accelStream = nullptr;
        decelStream = nullptr;

        FPSLimit = 500.0;

        gStream = 0;
        freqDelta = 0.0f;
        freqDelta2 = 0.0f;
        freqOldDelta = 0.0f;
        freqCurrent = -1.0f;
        freqMin = 0.0f;
        freqMax = 0.0f;
        accelVol = 1.0f;
        decelVol = 0.0f;

        rateVol = 1.0f;
        rateOldVol = rateVol;
        rateVolCurve = 10.0f;
        rateMinVol = 0.3f;
        rateEqAmount = 0.0f;
        rateEqOldAmount = 0.0f;
        rateEqCurve = 1.0f;
        rateRPMTarget = 1500.0f;
        rateOfChange = 0.0f;

        rateAccelXFadeRPMRatio = 2.0f; // 1/2 of the RPM range will be used as the crossfade range
        rateAccelXFadeRPMRange = 0.0f;
        rateDecelXFadeRPMRatio = 2.0f; // 1/2 of the RPM range will be used as the crossfade range
        rateDecelXFadeRPMRange = 0.0f;
        curRateRPMXFadeTarget = 0.0f;

        rateAccelEqXFadeRPMRatio = 4.0f; // 1/4 of the RPM range will be used as the crossfade range
        rateAccelEqXFadeRPMRange = 0.0f;
        rateDecelEqXFadeRPMRatio = 4.0f; // 1/4 of the RPM range will be used as the crossfade range
        rateDecelEqXFadeRPMRange = 0.0f;
        curRateEqRPMXFadeTarget = 0.0f;

        accelGlobalVol = 1.0f;
        decelGlobalVol = 1.0f;

        AccelXFadeRPMRatio = 50.0f; // 1/50 of the RPM range will be used as the crossfade range
        AccelXFadeRPMRange = 0.0f;
        DecelXFadeRPMRatio = 50.0f; // 1/50 of the RPM range will be used as the crossfade range
        DecelXFadeRPMRange = 0.0f;
        XFadeDelta = 0.0f;
        curXFadeRange = 0.0f;
        bAccelDirection = true;
        bOldDirection = true;
        bLoaded = false;
        
        hsIdle = 0;
        hsRedline = 0;
        hsReverseWhine = 0;
        hsForwardWhine = 0;
        hsIdleWhine = 0;

        chIdle = 0;
        chRedline = 0;
        chReverseWhine = 0;
        chForwardWhine = 0;
        chIdleWhine = 0;

        redlineStart = 0.0f;
        redlineVol = 0.0f;
        redlineGlobalVol = 1.0f;
        redlinePitch = 1.0f;
        redlineSampleRate = 0.0f;

        idleEnd = 0.0f;
        idleVol = 0.0f;
        idleGlobalVol = 1.0f;
        idlePitch = 1.0f;
        idleSampleRate = 0.0f;

        speedCurrent = 0.0f;
        speedDelta = 0.0f;
        speedOldDelta = 0.0f;

        reverseWhineFadeRange = 36.0f;
        reverseWhineMinVol = 0.5f;
        reverseWhineVol = 1.0f;
        reverseWhineGlobalVol = 1.0f;
        reverseWhineSampleRate = 0.0f;
        bReverseWhineEnable = false;

        forwardWhineFadeRange = 3.6f;
        forwardWhineMinVol = 0.1f;
        forwardWhineVol = 1.0f;
        forwardWhineGlobalVol = 1.0f;
        forwardWhineSampleRate = 0.0f;
        bForwardWhineEnable = false;

        idleWhineFadeRange = 7.2f;
        idleWhineMinVol = 0.0f;
        idleWhineVol = 1.0f;
        idleWhineGlobalVol = 1.0f;
        idleWhineSampleRate = 0.0f;
        bIdleWhineEnable = false;       

        bCurrentlyShifting = false;

        fxACL = 0;
        fxDCL = 0;

        eqMinGainACL = -15.0f;
        eqMinGainDCL = -15.0f;

        eqMaxGainACL = 0.0f;
        eqMaxGainDCL = 0.0f;

        eqFrequencyACL = 2756.25f;
        eqFrequencyDCL = 2756.25f;

        eqBandwidthACL = 36.0f;
        eqBandwidthDCL = 36.0f;

        eqACL = { 0 };
        eqACL.fBandwidth = eqBandwidthACL;
        eqACL.fCenter = eqFrequencyACL;
        eqACL.fGain = eqMaxGainACL;

        eqDCL = { 0 };
        eqDCL.fBandwidth = eqBandwidthDCL;
        eqDCL.fCenter = eqFrequencyDCL;
        eqDCL.fGain = eqMaxGainDCL;

        sampleRate = 0;
    }

    ~BASSGinsuMultiStream()
    {
        ReleaseEverything();
    }

};

class BASSGinsuPlayer
{
private:
    std::vector<BASSGinsuStream*> ginStreams;
    std::vector<BASSGinsuMultiStream*> ginMultiStreams;
    float volume;
    DWORD sampleRate;

public:
    BASSGinsuStream* CreateStream(std::filesystem::path ginPath)
    {
        BASSGinsuStream* stream = new BASSGinsuStream();
        if (stream == nullptr)
            return nullptr;

        if (!stream->Load(ginPath))
        {
            delete stream;
            return nullptr;
        }

        // initialize with global volume
        float v = stream->GetVolume();
        v = v * volume;
        stream->SetVolume(v);

        ginStreams.push_back(stream);
        return stream;
    }

    BASSGinsuMultiStream* CreateStream(std::filesystem::path accelGinPath, std::filesystem::path decelGinPath, std::filesystem::path redlinePath, std::filesystem::path idlePath, std::filesystem::path reverseWhinePath, std::filesystem::path forwardWhinePath, std::filesystem::path idleWhinePath)
    {
        BASSGinsuMultiStream* stream = new BASSGinsuMultiStream();
        if (stream == nullptr)
            return nullptr;

        if (!stream->Load(sampleRate, accelGinPath, decelGinPath, redlinePath, idlePath, reverseWhinePath, forwardWhinePath, idleWhinePath))
        {
            delete stream;
            return nullptr;
        }

        // initialize with global volume
        float v = stream->GetVolume();
        v = v * volume;
        stream->SetVolume(v);

        ginMultiStreams.push_back(stream);
        return stream;
    }

    void ReleaseStream(BASSGinsuStream* stream)
    {
        if (stream)
        {
            stream->Stop();

            auto it = std::find(ginStreams.begin(), ginStreams.end(), stream);
            if (it != ginStreams.end()) {
                ginStreams.erase(it);
            }

            delete stream;
        }
    }

    void ReleaseStream(BASSGinsuMultiStream* stream)
    {
        if (stream)
        {
            stream->Stop();

            auto it = std::find(ginMultiStreams.begin(), ginMultiStreams.end(), stream);
            if (it != ginMultiStreams.end()) {
                ginMultiStreams.erase(it);
            }

            delete stream;
        }
    }

    // Globally pause all streams.
    void Pause()
    {
        for (BASSGinsuStream* stream : ginStreams)
        {
            stream->Pause();
        }

        for (BASSGinsuMultiStream* stream : ginMultiStreams)
        {
            stream->Pause();
        }
    }

    // Globally play all streams.
    void Play()
    {
        for (BASSGinsuStream* stream : ginStreams)
        {
            stream->Play(FALSE);
        }
        for (BASSGinsuMultiStream* stream : ginMultiStreams)
        {
            stream->Play(FALSE);
        }
    }

    // Globally stop all streams.
    void Stop()
    {
        for (BASSGinsuStream* stream : ginStreams)
        {
            stream->Stop();
        }
        for (BASSGinsuMultiStream* stream : ginMultiStreams)
        {
            stream->Stop();
        }
    }

    // Sets global volume for all streams.
    void SetVolume(float inVol)
    {
        volume = inVol;

        for (BASSGinsuStream* stream : ginStreams)
        {
            float v = stream->GetVolume();
            v = v * volume;
            stream->SetVolume(v);
        }

        for (BASSGinsuMultiStream* stream : ginMultiStreams)
        {
            float v = stream->GetVolume();
            v = v * volume;
            stream->SetVolume(v);
        }
    }

    float GetVolume()
    {
        return volume;
    }

    BASSGinsuPlayer(int device, DWORD freq, DWORD flags, HWND win)
    {
        // BASS Init
        if (HIWORD(BASS_GetVersion()) != BASSVERSION) 
            throw std::runtime_error("An incorrect version of BASS.DLL was loaded.");
        if (!BASS_Init(device, freq, flags, win, NULL))
            throw std::runtime_error("BASS lib couldn't be initialized (BASS_Init)");
        
        BASS_SetConfig(BASS_CONFIG_MIXER_BUFFER, 0);
        BASS_SetConfig(BASS_CONFIG_MIXER_POSEX, 0);
        volume = 1.0f;
        sampleRate = freq;
    }

    ~BASSGinsuPlayer()
    {
        for (BASSGinsuStream* stream : ginStreams) 
        {
            delete stream;
        }
        for (BASSGinsuMultiStream* stream : ginMultiStreams)
        {
            delete stream;
        }

        BASS_Free();
    }
};

#endif