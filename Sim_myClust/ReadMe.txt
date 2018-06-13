## Synopsis
This is a Matlab demo codes of functional connectivity analysis based on the time-varying phase synchronization with a data-driven manner.
The basis of these programs are discribed described in the following paper (Yokoyama et al, 2018, https://www.frontiersin.org/articles/10.3389/fnhum.2018.00259/abstract).

See the details of the algrorithms in the programs for the description of the supplementary materials in above paper (correctly published
at this link: ghttps://www.frontiersin.org/articles/10.3389/fnhum.2018.
00259/full#supplementary-material.h).

## Source code
All programs are provided with Matlab format (*.m, or *.mexw32/64):

Main program:
	[dtw_c.c, dtw_c.mexw32, dtw_c.mexw64]
	  - These files are c-mex program for Dynamic time warping algorithms
   	 	See the following links for the details:
	            https://jp.mathworks.com/matlabcentral/fileexchange/43156-dynamic-time-warping--dtw-
@		@@https://wangquan.me/

	[Sim_myClust_demo.m] 
	  - This file is main matlab code of this simulation

Other functions, which are contained in the directory "/my_funtions":
	[cal_wpli.m]
	  - This is a function to calculate the phase-lag index (Vink, et, al. 2011)
	[figure_setting.m, figure_save.m]
	  - Configuration funtion for figure plot and save
        [Kuramoto_sim_v4.m]
          - This is a function to solve the Kuramoto model with "ode113" function. 
        [preprocessing.m]
          - This is a function for preprocessing of EEGs (including baseline correction, filtering, and z-score)
        [z_score.m]
          - This is a funtion to calculate the z-scored value

Requirements: 
 - EEGLAB toolbox;
   https://sccn.ucsd.edu/eeglab/index.php
   Download this toolbox, and add the path of saved directory of EEGLAB toolbox before use above scripts
   >> addpath('-saved directory of EEGLAB toolbox-')

## References

Yokoyama H, Nambu I, Izawa J and Wada Y. (2018). Alpha Phase Synchronization of Parietal Areas Reflects Switch-Specific Activity During Mental Rotation: an EEG study. Frontiers in Human Neuroscience,  doi: 10.3389/fnhum.2018.00259.

## License
Copyright: ? 2018 Yokoyama, Nambu, Izawa and Wada.
