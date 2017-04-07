
bids_dir='/Users/nickb/Analysis/fmrib/cvr_study_release/';

% load bids data structure
load([bids_dir 'derivatives/bids.mat']);

%add dependent paths
addpath([bids(1).dir 'derivatives/code/errorbarxy/']);

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

%produce figure s1
figure_plot_toronto_coherence(bids);

%produce figure s2
figure_plot_mean_etco2(bids);

%produce figure s3
figure_compare_tor_sin_ev(bids);

