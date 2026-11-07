# Intersection-point-height measure of human standing balance: computation and modeling
Supplementary material for "The intersection-point-height measure of human standing balance" (submitted)

## Quick Start Guide:
1. Ensure that you have MATLAB R2024a or later installed, along with the following toolboxes:
    > Control System Toolbox, Signal Processing Toolbox, Statistics and Machine Learning Toolbox.
2. Download standing data from Santos and Duarte (2016) and place the data files in the folder "Santos2016Data".
3. Run any of the scripts in the root folder to generate figures in the manuscript.

## Organization of the repository:

    \Forceplate_noise\ contains force plate measurement noise data

    \model_DIP\        contains functions related to the double-inverted-pendulum
                       balance model

    \Results\          contains analysis results that can be plotted
                       instead of computing the results from raw data
                       which can take some time

    \Santos2016Data\   download data from Santos and Duarte (2016) and add them
                       to this folder
                      
    \zIP\              contains functions for computing zIP from data 
                       and from models

## Acknowledgements:
Dr. Jongwoo Lee developed much of the original modeling code for the work published by Shiozawa et al. (2021).

## Patent Notice:
This software may be covered by one or more patents owned by Wisconsin Alumni Research Foundation.  No right or license to such patent(s) is granted.

## References:

Santos DA, Duarte M. A public data set of human balance evaluations. PeerJ 4: e2648, 2016. doi: 10.7717/peerj.2648.
  
Shiozawa K, Lee J, Russo M, Sternad D, Hogan N. Frequency-dependent force direction elucidates neural control of balance. J Neuroeng Rehabil 18: 145, 2021. doi: 10.1186/s12984-021-00907-2.

