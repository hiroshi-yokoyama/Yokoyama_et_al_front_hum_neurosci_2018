* Source code

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