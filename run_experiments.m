%% RUN_EXPERIMENTS  Controlled noise experiments (TEST data) + real noisy record.
%  Produces results/experiment_results.mat and results/tables/*.csv
%
%  Experiment design (see src/data_split.m for the exact train/test split)
%   - Ground truth : MIT-BIH records 100 and 103 (MLII), 60 s each, baseline removed
%   - Noise        : real NSTDB recordings: muscle artifact (MA), baseline
%                    wander (BW) and MA+BW (equal noise-power allocation)
%   - Noise level  : input SNR -5 ... 20 dB (PhysioNet nst definition)
%   - Monte Carlo  : 5 noise segments per condition and record (10 runs/point)
%   - Methods      : band-pass, wavelet, EKF, EKS (same baseline-removed input)


here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'));
if ~exist(fullfile(here, 'data', 'ecg_dataset.mat'), 'file'), prepare_dataset; end
ds = load(fullfile(here, 'data', 'ecg_dataset.mat'));
sp = data_split(ds);  fs = sp.fs;  L = sp.L;
outDir = fullfile(here, 'results');
if ~exist(fullfile(outDir, 'tables'), 'dir'), mkdir(fullfile(outDir, 'tables')); end

% EKF hyper-parameters: from tune_parameters.m if it has been run
tp = fullfile(outDir, 'tuned_params.mat');
if exist(tp, 'file'), load(tp, 'opt'); fprintf('Using tuned parameters from %s\n', tp);
else, opt = struct; fprintf('tuned_params.mat not found -> using defaults in ekf_ecg.m\n'); end

%% ---------- 1. Noise power (PhysioNet nst rules) ----------
[~, N_ma] = nst_power([], [], [], fs, ds.noise_ma);
[~, N_bw] = nst_power([], [], [], fs, ds.noise_bw);

%% ---------- 2. Clean TEST segments ----------
nRec = numel(sp.test.rec);
for q = 1:nRec
    [rec(q).xc, rec(q).pkRef, Sq] = get_clean_segment(ds, sp.test.rec(q), sp.test.start(q), L);
    rec(q).P = struct('S', Sq, 'N_ma', N_ma, 'N_bw', N_bw);
    rec(q).id = sp.test.rec(q);
    fprintf('Record %d: nst S = %.4f mV^2, %d beats in test segment\n', ...
            rec(q).id, Sq, numel(rec(q).pkRef));
end

%% ---------- 3. Demo case: record 100, MA+BW noise at 5 dB ----------
demo.type = 'ma+bw';  demo.snr = 5;  demo.offset = sp.test.noiseOffsets(1);
demo.xc = rec(1).xc;
demo.v  = make_noise(ds, demo.type, demo.snr, L, demo.offset, rec(1).P);
demo.y  = demo.xc + demo.v;
demo.r  = ekf_pipeline(demo.y, fs, opt);
demo.bp = denoise_bandpass(demo.r.xb, fs);
demo.wt = denoise_wavelet(demo.r.xb);
est = {demo.bp, demo.wt, demo.r.xEKF, demo.r.xEKS};
for j = 1:4, demo.m(j) = eval_metrics(demo.xc, demo.y, est{j}, fs); end
demo.m0 = eval_metrics(demo.xc, demo.y, demo.r.xb, fs);   % baseline-removed input
fprintf('Demo (%s, %d dB): output SNR  BP %.2f | WT %.2f | EKF %.2f | EKS %.2f dB\n', ...
        demo.type, demo.snr, [demo.m.snrOut]);

%% ---------- 4. SNR sweep (Monte Carlo) ----------
types   = {'ma', 'bw', 'ma+bw'};
snrs    = -5:2.5:20;
offsets = sp.test.noiseOffsets;  nTrial = numel(offsets);
methods = {'bandpass', 'wavelet', 'ekf', 'eks'};
sz = [nRec, numel(types), numel(snrs), nTrial];
S = struct('snrIn', nan(sz), 'snrBL', nan(sz), 'snrOut', nan([sz 4]), 'snrImp', nan([sz 4]), ...
           'rmse', nan([sz 4]), 'prd', nan([sz 4]), 'corr', nan([sz 4]), ...
           'sens', nan(sz), 'ppv', nan(sz));
tic;
for q = 1:nRec
  for a = 1:numel(types)
    for b = 1:numel(snrs)
      for t = 1:nTrial
        y = rec(q).xc + make_noise(ds, types{a}, snrs(b), L, offsets(t), rec(q).P);
        r = ekf_pipeline(y, fs, opt);
        est = {denoise_bandpass(r.xb, fs), denoise_wavelet(r.xb), r.xEKF, r.xEKS};
        for j = 1:4
            m = eval_metrics(rec(q).xc, y, est{j}, fs);
            S.snrOut(q,a,b,t,j) = m.snrOut;  S.snrImp(q,a,b,t,j) = m.snrImp;
            S.rmse(q,a,b,t,j)   = m.rmse;    S.prd(q,a,b,t,j)    = m.prd;
            S.corr(q,a,b,t,j)   = m.corr;
        end
        S.snrIn(q,a,b,t) = m.snrIn;
        m0 = eval_metrics(rec(q).xc, y, r.xb, fs);  S.snrBL(q,a,b,t) = m0.snrOut;
        [S.sens(q,a,b,t), S.ppv(q,a,b,t)] = peak_detection_stats(r.pk, rec(q).pkRef, fs);
      end
      fprintf('rec %d %-6s %5.1f dB  out-SNR: BP %5.2f WT %5.2f EKF %5.2f EKS %5.2f  (Se %.2f PPV %.2f)\n', ...
          rec(q).id, types{a}, snrs(b), squeeze(mean(S.snrOut(q,a,b,:,:), 4)), ...
          mean(S.sens(q,a,b,:)), mean(S.ppv(q,a,b,:)));
    end
  end
end
fprintf('Sweep finished in %.1f s\n', toc);

%% ---------- 5. Tables (CSV): pooled over both test records ----------
pool = @(X, a, b, j) X(:,a,b,:,j);
rows = {};
for a = 1:numel(types)
    for b = 1:numel(snrs)
        for j = 1:4
            v = pool(S.snrOut, a, b, j);  w = pool(S.snrImp, a, b, j);
            rows(end+1,:) = {types{a}, snrs(b), methods{j}, ...
                mean(S.snrIn(:,a,b,:), 'all'), mean(v(:)), std(v(:)), mean(w(:)), ...
                mean(pool(S.rmse,a,b,j), 'all'), mean(pool(S.prd,a,b,j), 'all'), ...
                mean(pool(S.corr,a,b,j), 'all'), mean(S.sens(:,a,b,:), 'all'), ...
                mean(S.ppv(:,a,b,:), 'all')}; %#ok<SAGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'noise', 'snr_nst_dB', 'method', 'snr_in_dB', ...
    'snr_out_dB', 'snr_out_std_dB', 'snr_imp_dB', 'rmse_mV', 'prd_pct', 'corr', ...
    'rpeak_sens', 'rpeak_ppv'});
writetable(T, fullfile(outDir, 'tables', 'snr_sweep_results.csv'));

%% ---------- 6. Real recording: MIT-BIH record 105 (no ground truth) ----------
% (a) 20:00-21:00 : the noisiest part of the record (dense noise annotations)
% (b) 5:00-6:00   : cleaner part containing premature ventricular beats (PVCs),
%                   used to show a limitation of the model-based filter
segStart = [20*60, 5*60];
for q = 1:2
    j0 = segStart(q)*fs;
    rr.y  = ds.ecg105(j0+1:j0+L);
    rr.r  = ekf_pipeline(rr.y, fs, opt);
    rr.bp = denoise_bandpass(rr.r.xb, fs);
    rr.wt = denoise_wavelet(rr.r.xb);
    rr.t0 = j0/fs;
    a = ds.ann105_samp;  k = a > j0 & a <= j0 + L;
    rr.annSamp = a(k) - j0;  rr.annType = ds.ann105_type(k);
    real(q) = rr; %#ok<SAGROW>
end

save(fullfile(outDir, 'experiment_results.mat'), 'rec', 'demo', 'S', 'types', 'snrs', ...
     'nTrial', 'offsets', 'methods', 'real', 'fs', 'L', 'opt', 'sp');
fprintf('Saved results/experiment_results.mat\n');
