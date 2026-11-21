%% TUNE_PARAMETERS  Choose the EKF process-noise hyper-parameters (TRAIN data only).
%  Grid search that maximises the mean EKS output SNR over 5 training cases
%  (record 100, 1:00-2:00, fixed noise segment).  Saves results/tuned_params.mat,
%  which run_experiments.m uses. 

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'));
if ~exist(fullfile(here, 'data', 'ecg_dataset.mat'), 'file'), prepare_dataset; end
ds = load(fullfile(here, 'data', 'ecg_dataset.mat'));
sp = data_split(ds);  fs = sp.fs;  L = sp.L;
if ~exist(fullfile(here, 'results'), 'dir'), mkdir(fullfile(here, 'results')); end

[xc, ~, P.S] = get_clean_segment(ds, sp.train.rec, sp.train.start, L);
[~, P.N_ma] = nst_power([], [], [], fs, ds.noise_ma);
[~, P.N_bw] = nst_power([], [], [], fs, ds.noise_bw);

cases = {'ma', 5; 'ma', 10; 'bw', 10; 'ma+bw', 5; 'ma', 20};   % training conditions
Y = cell(size(cases,1), 1);
for c = 1:size(cases,1)
    Y{c} = xc + make_noise(ds, cases{c,1}, cases{c,2}, L, sp.train.noiseOffset, P);
end

qgrid.qa   = [0.05 0.1 0.2];       % std of Gaussian amplitudes (fraction of alpha)
qgrid.qth  = [0.01 0.02 0.05];     % std of Gaussian centres (rad)
qgrid.qeta = [0.001 0.002 0.005 0.01];  % std of additive model noise (fraction of max alpha)
best.score = -inf;
res = [];
for qa = qgrid.qa
    for qth = qgrid.qth
        for qeta = qgrid.qeta
            o = struct('qa', qa, 'qth', qth, 'qeta', qeta, 'qb', 0.05, 'qw', 0.05);
            v = zeros(1, size(cases,1));
            for c = 1:size(cases,1)
                r = ekf_pipeline(Y{c}, fs, o);
                m = eval_metrics(xc, Y{c}, r.xEKS, fs);
                v(c) = m.snrOut;
            end
            res(end+1,:) = [qa qth qeta mean(v)]; %#ok<SAGROW>
            fprintf('qa=%.3f  qth=%.3f  qeta=%.4f  ->  mean EKS SNR %.2f dB\n', qa, qth, qeta, mean(v));
            if mean(v) > best.score, best.score = mean(v); best.opt = o; end
        end
    end
end
opt = best.opt;
save(fullfile(here, 'results', 'tuned_params.mat'), 'opt', 'res', 'qgrid', 'cases');
fprintf('\nBest: qa=%.3f qth=%.3f qeta=%.4f (mean EKS SNR %.2f dB)\n', ...
        opt.qa, opt.qth, opt.qeta, best.score);
