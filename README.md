# HybTrack_for_TS
We modified the HybTrack MATLAB script for transcription site detection
The original code is from http://rdcu.be/EnAp Lee, B. H. and Park, H. Y. (2018). HybTrack: A Hybrid Single Particle Tracking Software Using Manual and Automatic Detection of Dim Signals. Scientific Reports, 8, 212.
Now modified particle detection algorithm sorts the pixels' intensity from the croped image as 1d array.
This sorted intensities are fitted into inverse error function. The algorithm calculates ratio between rms values of fitting errors from left and right side of sorted intensities.
If the ratio is greater than threshold, it is discriminated that the croped image has transcription site.

Particle selection from the first frame was changed as maximum projected image
