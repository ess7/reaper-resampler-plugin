# reaper-resampler-plugin

Replaces Reaper's resampler by hooking [`WDL_Resampler::ResampleOut`](https://github.com/justinfrankel/WDL/blob/master/WDL/resample.cpp). Quality is similar to SoX VHQ.

[Sine sweep spectrum comparisons](https://imgur.com/a/A3rB3iT)

### Installation

Copy `reaper_resample32.dll` (for 32-bit) or `reaper_resample64.dll` (for 64-bit) to Reaper's `Plugins` folder. You should see a "Resampler loaded" message box in the background. Otherwise, create an issue with the Reaper version/OS you are using.

#### Tested versions (32/64-bit):
* 4.77 - 4.78
* 5.941 - 5.95
