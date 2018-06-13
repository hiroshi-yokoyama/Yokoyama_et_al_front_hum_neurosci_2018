* Data set <Sim_myClust>

Simulation data are provided in Matlab format (*.mat):

[rnd_seed.mat]
  - This is the setting file of random numner generator, which are applied to generate the simulated EEG signals with the random noise in main script "Sim_myClust.m"
    If you run "Sim_myClust.m" with the settings of random generator in this mat-file, 
    the simulation results, which were mentioned in the Supplementary Materials, can be reconstructed in any time. 

[Sim_myClust_dataset.mat 
  - This file is contained all analysed data for the simulation, which were mentioned in the Supplementary Materials.

All data sed with Matlab fromat containing variables :
    + data        ;  simulated EEG data, 3-D matrix [ch x time x trial]
    + t           ;  [sample x 1] colmun vector indicating the time periods
    
    + theta1      ;  phase angle for each oscillator in oscillator group 1 [Nosc x time sample]
    + theta2      ;  phase angle for each oscillator in oscillator group 2 [Nosc x time sample]
    + signal1     ;  time-domain signals for each oscillator in oscillator group 1 [Nosc x time sample]
    + signal2     ;  time-domain signals for each oscillator in oscillator group 2 [Nosc x time sample]
    + Tout        ;  out put variabe of "Kuramoto_sim_v4()" function 
                     [sample x 1] colmun vector indicating the time periods
                     see the detail in descriptions of "Kuramoto_sim_v4()" function
    + rnd_seed1   ;  the settings of the random number generator for osicllator group 1
    + rnd_seed2   ;  the settings of the random number generator for osicllator group 2
    + rnd_seed3   ;  the settings of the random number generator for random noise in each trial
    
    
    + alltf       ;  time frequency of simulated data. [freqs x time x trial x ch] 
                     the data structure is in accordance with the EEGLAB functions  (See the details in descriptions of timefreq() function in EEGLAB)
    + freqs       ;  vector of frequency band, size [1 x frequeancy]
                     output variable of "timefreq()" function in EEGLAB
    + t_tf        ;  vector of time periods, size [1 x time]
                     output variable of "timefreq()" function in EEGLAB
    
    + wpli        ;  time-varying weighted phase lag indeces (wPLIs) each pair of simulated EEG signals (z-scored)
    + clust       ;  struture of clustering result
                     - clust.DTW
		       the data structure is in accordance with the MATLAB functions  (See the details in descriptions of cluster() function in MATLAB)
    + eva_DTW     ;  struture of evaluation result of optimal cluster number
                       the data structure is in accordance with the MATLAB functions  (See the details in descriptions of evalclusters() function in MATLAB)