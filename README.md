# Ratiomtric-4Pi
 Ratimotric-4Pi is a graphics processing unit (GPU) based global fitting algorithm for 4Pi-SMLM with flexible PSF modeling and parameter sharing, to extract maximum information from 4Pi single molecule data and achieved both good color separation and optimal 3D resolution. By partially linking the photon parameters between channels with interference difference of π during global fitting of the multi-channel 4Pi single molecule data, we showed on simulated data that the loss of the localization precision is minimal compared with the theoretical minimum uncertainty, the Cramer-Rao lower bound (CRLB). The fitting speeds achieve ~5845 fits/s on a standard GPU (NVIDIA RTX3090) for regions of interest (ROI) with a size of 13×13 pixels.
 
 <img src="https://github.com/Li-Lab-SUSTech/Ratiometric-4Pi/blob/main/Figure/Fig_1_Schematic of 4Pi-SMLM.png" width = 60%  alt="workflow overview" align=center />

# Requirements
Matlab R2019a or newer  

The GPU fitter requires:
  
  - Microsoft Windows 7 or newer, 64-bit
  - CUDA capable graphics card, minimum Compute Capability 6.1
  - CUDA 11.3 compatible graphics driver (for GeForce products 471.41 or later)

The CPU version runs on macOS and Microsoft Windows 7 or newer, 64-bit

 # How to run
 Examples code are avalible in file **multicolour_4Pi_random_shared.m, multicolour_CRLB_3PATH.m, multicolour_CRLB_4PATH.m, multicolour_CRLB_Comparison_YD.m**. 
 
# Contact
For any questions / comments about this software, please contact [Li Lab](https://faculty.sustech.edu.cn/liym2019/en/).

# Copyright 
Copyright (c) 2021 Li Lab, Southern University of Science and Technology, Shenzhen.
