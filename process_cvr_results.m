%add dependent paths
addpath ~/Documents/MATLAB/json4mat
addpath /usr/local/fsl/etc/matlab

bids_dir='/Users/nickb/Analysis/fmrib/cvr_study/';

% load bids data structure
load([bids_dir 'derivatives/bids.mat']);

close all;

%produce table 1
T1=tabulate_respdata(bids);

%produce table 2
T2=tabulate_stats(bids);

%produce figure 2
figure_cvrmaps(bids,1);

%produce figure 3
figure_compare_tor_sin(bids);

%produce figure 4
figure_sincvrmaps(bids,1);

%produce figure 5
figure_compare_sin_sin(bids);

