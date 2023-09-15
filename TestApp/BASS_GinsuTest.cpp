// BASSGinsu Test App code
// Very jank.
//

#define NOMINMAX
#include <iostream>
#include <thread>
#include <windows.h>
#include "../BASSGinsu.hpp"
#include <Xinput.h>

#pragma comment(lib, "XInput.lib")

#define RPMCHANGE 3.0f
#define SPEEDCHANGE 0.05f

float testFreq;
float testFreqMax;
float testFreqMin;

float testSpeed = 0.0f;
float testSpeedMax = 300.0f;

constexpr bool bIsFloat = false;

BASSGinsuMultiStream* ginMultiStream;

HSTREAM testStream;

float GetRightTriggerValue(DWORD controllerIndex) {
    XINPUT_STATE state;
    ZeroMemory(&state, sizeof(XINPUT_STATE));

    // Get the state of the controller
    if (XInputGetState(controllerIndex, &state) == ERROR_SUCCESS) {
        // Normalize the trigger value to a float between 0 and 1.0
        float rightTrigger = static_cast<float>(state.Gamepad.bRightTrigger) / 255.0f;
        return rightTrigger;
    }
    else {
        // Controller not connected or an error occurred
        return 0.0f;
    }
}

void HandleButtonPress()
{
    float throttleRT = GetRightTriggerValue(0);
    float adder = cus_lerp(0.0f, RPMCHANGE, throttleRT);

    if (adder > 0)
        testFreq += adder;
    else if (GetAsyncKeyState(VK_LCONTROL) >> 15)
    {
        testFreq += RPMCHANGE;
        testSpeed += SPEEDCHANGE;
    }
    else
    {
        testFreq -= RPMCHANGE;
        testSpeed -= SPEEDCHANGE;
    }


    float adder2 = cus_lerp(0.0f, SPEEDCHANGE, throttleRT);
    if (adder2 > 0)
        testSpeed += adder2;
    else if (GetAsyncKeyState(VK_LCONTROL) >> 15)
    {
        testSpeed += SPEEDCHANGE;
    }
    else
    {
        testSpeed -= SPEEDCHANGE;
    }

    testFreq = std::clamp(testFreq, testFreqMin, testFreqMax);
    testSpeed = std::clamp(testSpeed, 0.0f, testSpeedMax);
}


float deltaChange = RPMCHANGE;
float bounceFloat(float& currentValue, float minValue, float maxValue) {
    // Check if the current value is within the bounds
    if (currentValue >= minValue && currentValue <= maxValue) {
        // Update the current value based on the delta
        currentValue += deltaChange;

        // Check if the new value is within bounds
        if (currentValue < minValue) {
            currentValue = minValue;
            deltaChange = -deltaChange;  // Reverse the direction
        }
        else if (currentValue >= maxValue) {
            currentValue = maxValue;
            deltaChange = -deltaChange;  // Reverse the direction
        }
    }
    else {
        // Reset the current value to the nearest bound
        if (currentValue < minValue) {
            currentValue = minValue;
        }
        else {
            currentValue = maxValue;
        }
        deltaChange = -deltaChange;  // Reverse the direction
    }

    return currentValue;
}

void BounceRPM()
{
    bounceFloat(testFreq, testFreqMin, testFreqMax);
}

void inputloop()
{
    while (1)
    {
        HandleButtonPress();
        //BounceRPM();
        Sleep(1);
    }
}

DWORD CALLBACK StreamProc(HSTREAM handle, void* buffer, DWORD length, void* user)
{
    ginMultiStream->GetData((int16_t*)buffer, length);
    return length;
}

void thyLoop()
{
    ginMultiStream->SetFrequency(testFreq);
    ginMultiStream->SetSpeed(testSpeed / 3.6f);
    Sleep(1);
}

int main(int argc, char** argv)
{
    std::cout << "BASS Ginsu Synthesizer Test App\nHold LCONTROL for throttle\nOr use RT on a XInput controller as throttle\n";

    BASSGinsuPlayer ginPlayer(-1, 44100, 0, NULL);

    ginMultiStream = ginPlayer.CreateStream(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], bIsFloat);

    // test values from NFSMW
    testFreqMin = 1500.0f;
    testFreqMax = 7784.0f;

    ginMultiStream->SetMinFrequency(testFreqMin);
    ginMultiStream->SetMaxFrequency(testFreqMax);
    testFreq = testFreqMin;

    int d = std::stoi(argv[8]);
    switch (d)
    {
    case 2:
        ginMultiStream->SetIdleWhineEnable(true);
        ginMultiStream->SetForwardWhineEnable(true);
        break;
    case 1:
        ginMultiStream->SetReverseWhineEnable(true);
        break;
    default:
        break;
    }

    

    if (bIsFloat)
        testStream = BASS_StreamCreate(ginMultiStream->GetSampleRate(), 1, BASS_SAMPLE_FLOAT, StreamProc, 0);
    else
        testStream = BASS_StreamCreate(ginMultiStream->GetSampleRate(), 1, 0, StreamProc, 0);

    BASS_ChannelSetAttribute(testStream, BASS_ATTRIB_BUFFER, 0);
    BASS_ChannelPlay(testStream, FALSE);

    std::thread threed(inputloop);
    threed.detach();


    while (1)
    {
        thyLoop();
    }

    getchar();

    return 0;
}
