%% PREPARE_DATASET  Convert raw PhysioNet WFDB files into one MATLAB dataset.
%  Run once.  Creates data/ecg_dataset.mat, which main.m uses.
%
%  Sources (PhysioNet, Open Data Commons Attribution License):
%   - MIT-BIH Arrhythmia Database   : records 100 & 103 (clean ECG), 105 (noisy ECG)
%   - MIT-BIH Noise Stress Test DB   : 'ma' (muscle artifact), 'bw' (baseline wander)
%  All records are sampled at 360 Hz.

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'));
rawDir = fullfile(here, 'data', 'raw');

% ---------- clean ECG: record 100, lead MLII ----------
[s100, fs] = read_wfdb212(fullfile(rawDir, '100'));
[annSamp, annType] = read_wfdb_ann(fullfile(rawDir, '100.atr'));
ds.fs           = fs;
ds.ecg100       = s100(:,1);              % MLII lead, mV
ds.ann100_samp  = annSamp;
ds.ann100_type  = annType;                % 1 = normal beat (N)

% ---------- second clean ECG (different patient): record 103, lead MLII ----------
s103 = read_wfdb212(fullfile(rawDir, '103'));
[a103s, a103t] = read_wfdb_ann(fullfile(rawDir, '103.atr'));
ds.ecg103      = s103(:,1);
ds.ann103_samp = a103s;
ds.ann103_type = a103t;

% ---------- naturally noisy ECG: record 105, lead MLII ----------
s105 = read_wfdb212(fullfile(rawDir, '105'));
[a105s, a105t] = read_wfdb_ann(fullfile(rawDir, '105.atr'));
ds.ecg105      = s105(:,1);
ds.ann105_samp = a105s;
ds.ann105_type = a105t;

% ---------- real noise recordings (channel 1) ----------
ma = read_wfdb212(fullfile(rawDir, 'ma'));
bw = read_wfdb212(fullfile(rawDir, 'bw'));
ds.noise_ma = ma(:,1);
ds.noise_bw = bw(:,1);

ds.description = ['MIT-BIH records 100, 103 & 105 (MLII) + NSTDB ma/bw noise (ch 1), ' ...
                  'fs = 360 Hz, units mV.  Source: physionet.org'];

save(fullfile(here, 'data', 'ecg_dataset.mat'), '-struct', 'ds', '-v7');
fprintf('Saved data/ecg_dataset.mat  (%d samples per record, fs = %d Hz)\n', ...
        numel(ds.ecg100), fs);
