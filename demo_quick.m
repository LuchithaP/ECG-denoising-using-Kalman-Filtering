%% DEMO_QUICK  Interactive demo: denoise one noisy ECG and show the result.
%   Change the three settings below and press Run.

noiseType = 'ma+bw';   % 'ma' (muscle), 'bw' (baseline wander) or 'ma+bw'
snrDb     = 5;         % input SNR in dB (PhysioNet nst definition)
record    = 100;       % clean ECG: 100 or 103

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'));
ds = load(fullfile(here, 'data', 'ecg_dataset.mat'));
sp = data_split(ds);  fs = ds.fs;  L = sp.L;

% ---- build the noisy test signal (ground truth known) ----
[xc, ~, S] = get_clean_segment(ds, record, 300*fs, L);
[~, Nma] = nst_power([], [], [], fs, ds.noise_ma);
[~, Nbw] = nst_power([], [], [], fs, ds.noise_bw);
P = struct('S', S, 'N_ma', Nma, 'N_bw', Nbw);
y = xc + make_noise(ds, noiseType, snrDb, L, sp.test.noiseOffsets(1), P);

% ---- estimate ----
tic; r = ekf_pipeline(y, fs); tRun = toc;      % LS fit + MLE variance + EKF + EKS
bp = denoise_bandpass(r.xb, fs);
wt = denoise_wavelet(r.xb);

% ---- metrics ----
names = {'Band-pass', 'Wavelet', 'EKF', 'EKS'};  est = {bp, wt, r.xEKF, r.xEKS};
fprintf('\n%s noise, input SNR %g dB, record %d  (pipeline %.2f s)\n', noiseType, snrDb, record, tRun);
fprintf('%-10s %10s %10s %8s\n', 'Method', 'SNRout(dB)', 'RMSE(mV)', 'corr');
for j = 1:4
    m = eval_metrics(xc, y, est{j}, fs);
    fprintf('%-10s %10.2f %10.3f %8.3f\n', names{j}, m.snrOut, m.rmse, m.corr);
end

% ---- plot ----
C = fig_style();  t = (0:L-1)'/fs;  k = t >= 20 & t < 25;
figure('Color', 'w', 'Position', [50 50 1600 800]);
subplot(2,1,1); plot(t(k), y(k), 'Color', C.noisy); hold on; plot(t(k), xc(k), 'k', 'LineWidth', 1.5);
title(sprintf('Noisy input (%s, %g dB)', noiseType, snrDb), 'FontSize', 16); legend('Noisy', 'Clean'); grid on;
subplot(2,1,2); plot(t(k), xc(k), 'k', 'LineWidth', 1.5); hold on;
plot(t(k), bp(k), 'Color', C.bp, 'LineWidth', 1.5); plot(t(k), r.xEKS(k), 'Color', C.eks, 'LineWidth', 2);
title('Estimates', 'FontSize', 16); legend('Clean', 'Band-pass', 'EKS'); xlabel('Time (s)'); grid on;
