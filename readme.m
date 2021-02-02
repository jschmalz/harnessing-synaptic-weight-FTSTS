% The following matlab files generate the weight profiles for Figures 3 - 14 
% in the paper titled "Controlling Synchronization of Spiking Neuronal Networks by Harnessing
% Synaptic Plasticity" by Joseph Schmalz and Gautam Kumar. A rastor plot of
% the spike times can be generated by uncommenting spike_E and spike_I in
% the main_ matlab file code. Please note that recording the spike times will increase the
% run time and memory requirement of the code.
%
% Each folder contains three matlab files. The file entitled main_... initiates
% the code. The file entitle ode_... contains the ODEs for the system. The
% file entitle pulsatile_... creates the FTSTS or CR pulse applied to the
% network.
%
% In folders where multiple figures are listed in the folder name, different sections of the
% code must be comment/uncomment and are designated by "TO MAKE ...". These can be found in BOTH the main_ and
% ode_ files so make sure to look in both. 