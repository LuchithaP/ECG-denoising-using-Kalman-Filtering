%% MAIN  ECG denoising with Least Squares + Extended Kalman Filter/Smoother
%  Run this file to reproduce every result.  Each step can also be run on its own.
%
%   Step 1  prepare_dataset   raw PhysioNet files  -> data/ecg_dataset.mat     
%   Step 2  tune_parameters   EKF hyper-parameters on TRAIN data             
%   Step 3  run_experiments   TEST: SNR sweep, demo case, record 105         
%   Step 4  make_figures      results/figures/*.png  (1920x1080)             
%   Step 5  make_videos       results/videos/*.mp4   (1920x1080, 30 fps)    
%   Step 6  demo_quick        one noisy ECG, plots on screen (for live demo) 
%
%  Requirements:  Signal Processing Toolbox,
%                 Wavelet Toolbox (only for the wavelet baseline).

clear; clc; close all;
here = fileparts(mfilename('fullpath'));
cd(here);  addpath(fullfile(here, 'src'));

prepare_dataset;      % Step 1
tune_parameters;      % Step 2  (optional: defaults in src/ekf_ecg.m are the tuned values)
run_experiments;      % Step 3
make_figures;         % Step 4
make_videos;          % Step 5
