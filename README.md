# reaper-resampler-plugin

Replaces Reaper's resampler by hooking [`WDL_Resampler::ResampleOut`](https://github.com/justinfrankel/WDL/blob/master/WDL/resample.cpp).
Optimized for conversion between common sample rates (44.1, 48, 88.2, 96, 176.4, 192).
Quality is similar to SoX VHQ.

[Sine sweep spectrum comparisons](https://imgur.com/a/A3rB3iT)

[CPU usage tests](https://github.com/ess7/reaper-resampler-plugin/wiki/CPU-usage-tests)

### Installation

Copy `reaper_resample32.dll` (for 32-bit) or `reaper_resample64.dll` (for 64-bit) to Reaper's `Plugins` folder. You should see a "Resampler loaded" message box in the background. Otherwise, create an issue with the Reaper version/OS you are using.

#### Tested versions (32/64-bit):
* 4.77 - 4.78
* 5.941 - 5.95
