%% MAKE_FIGURES  All result figures as 1920x1080 (16:9) PNGs 
%  Requires results/experiment_results.mat (run run_experiments.m first).
%  Output: results/figures/fig01 ... fig09

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'));
ds = load(fullfile(here, 'data', 'ecg_dataset.mat'));
R  = load(fullfile(here, 'results', 'experiment_results.mat'));
figDir = fullfile(here, 'results', 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end
C  = fig_style();
fs = R.fs;  d = R.demo;  t = (0:R.L-1)'/fs;
tw = @(a, b) t >= a & t < b;                 % time-window helper

%% ===== Fig 1: the problem - clean ECG corrupted by real noise =====
[~, f] = fig_style(true);
vMA = make_noise(ds, 'ma', d.snr + 10*log10(2), R.L, d.offset, R.rec(1).P);  % MA part
vBW = make_noise(ds, 'bw', d.snr + 10*log10(2), R.L, d.offset, R.rec(1).P);  % BW part
k = tw(20, 26);
tl = tiledlayout(f, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
sigs = {d.xc, vMA, vBW, d.y};
ttl = {'Clean ECG  (MIT-BIH record 100, lead MLII)', ...
       'Muscle artifact (real EMG noise, NSTDB "ma")', ...
       'Baseline wander (real motion/respiration noise, NSTDB "bw")', ...
       sprintf('Recorded signal = ECG + MA + BW   (input SNR = %d dB, PhysioNet nst definition)', d.snr)};
col = {C.clean, C.wt, C.ekf, C.noisy*0.7};
for i = 1:4
    ax = nexttile(tl); plot(ax, t(k), sigs{i}(k), 'Color', col{i}, 'LineWidth', 1.6);
    title(ax, ttl{i}, 'FontSize', 17); ylabel(ax, 'mV'); grid(ax, 'on'); xlim(ax, [20 26]);
    if i < 4, ax.XTickLabel = []; end
end
xlabel(ax, 'Time (s)');
grab_frame(f, fullfile(figDir, 'fig01_problem_noisy_ecg.png')); close(f);

%% ===== Fig 2: dynamical ECG model - phase and limit cycle =====
[~, f] = fig_style(true);
tl = tiledlayout(f, 2, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
k = tw(20, 23);  pk = d.r.pk(d.r.pk > find(k,1) & d.r.pk < find(k,1,'last'));
ax = nexttile(tl, 1);
plot(ax, t(k), d.xc(k), 'k', 'LineWidth', 2); hold(ax, 'on');
plot(ax, t(pk), d.xc(pk), 'v', 'MarkerSize', 12, 'MarkerFaceColor', C.eks, 'MarkerEdgeColor', C.eks);
title(ax, 'ECG with detected R peaks'); ylabel(ax, 'mV'); grid(ax, 'on');
ax = nexttile(tl, 3);
plot(ax, t(k), d.r.phi(k), 'Color', C.ekf, 'LineWidth', 2.2);
yticks(ax, [-pi 0 pi]); yticklabels(ax, {'-\pi', '0', '\pi'});
title(ax, 'Cardiac phase \theta(t): one 2\pi turn per beat');
xlabel(ax, 'Time (s)'); ylabel(ax, '\theta (rad)'); grid(ax, 'on');
ax = nexttile(tl, 2, [2 1]);
kk = tw(10, 30);
plot3(ax, cos(d.r.phi(kk)), sin(d.r.phi(kk)), d.xc(kk), 'Color', [C.ekf 0.5], 'LineWidth', 1.2); hold(ax, 'on');
th = linspace(-pi, pi, 400);
plot3(ax, cos(th), sin(th), gauss_ecg_model(d.r.p, th'), 'Color', C.eks, 'LineWidth', 3.5);
plot3(ax, cos(th), sin(th), 0*th, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
view(ax, -35, 22); grid(ax, 'on');
xlabel(ax, 'cos \theta'); ylabel(ax, 'sin \theta'); zlabel(ax, 'z (mV)');
title(ax, {'State space: the ECG is a point moving', 'around a circle (limit cycle)'});
legend(ax, {'20 beats (clean)', 'LS-fitted Gaussian model'}, 'Location', 'southoutside');
grab_frame(f, fullfile(figDir, 'fig02_dynamic_model.png')); close(f);

%% ===== Fig 3: Least-squares fit of the 5-Gaussian beat model =====
[~, f] = fig_style(true);
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl, 1, [2 2]); hold(ax, 'on');
sub = 1:7:R.L;
scatter(ax, d.r.phi(sub), d.r.xb(sub), 6, C.noisy, 'filled', 'MarkerFaceAlpha', 0.25);
th = linspace(-pi, pi, 600)';
n = 5; p = d.r.p; wn = {'P', 'Q', 'R', 'S', 'T'};
for i = 1:n
    gi = p(i)*exp(-wrap_pi(th - p(2*n+i)).^2/(2*p(n+i)^2));
    fill(ax, [th; flipud(th)], [gi; zeros(size(gi))], C.waves(i,:), ...
         'FaceAlpha', 0.35, 'EdgeColor', C.waves(i,:), 'LineWidth', 1.5);
    [~, j] = max(abs(gi));
    text(ax, th(j), gi(j) + 0.12*sign(gi(j)), wn{i}, 'FontSize', 22, 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'center', 'Color', C.waves(i,:)*0.8);
end
h1 = plot(ax, d.r.phc, d.r.mu, 'k', 'LineWidth', 3);
h2 = plot(ax, th, gauss_ecg_model(p, th), '--', 'Color', C.eks, 'LineWidth', 3);
xlim(ax, [-pi pi]); xticks(ax, -pi:pi/2:pi);
ylim(ax, [min(d.r.mu) - 0.35, max(d.r.mu) + 0.35]);
xticklabels(ax, {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
xlabel(ax, 'Cardiac phase \theta (rad)'); ylabel(ax, 'Amplitude (mV)'); grid(ax, 'on');
legend(ax, [h1 h2], {'Phase-averaged noisy beat', 'LS fit: sum of 5 Gaussians'}, 'Location', 'northwest');
title(ax, 'Least squares: fit z(\theta) = \Sigma_i \alpha_i exp(-(\theta-\theta_i)^2 / 2b_i^2) to the noisy data');
ax = nexttile(tl, 3);
semilogy(ax, 0:numel(d.r.lsHist.cost)-1, d.r.lsHist.cost, 'o-', 'Color', C.ekf, ...
         'MarkerFaceColor', C.ekf, 'LineWidth', 2);
xlabel(ax, 'Levenberg-Marquardt iteration'); ylabel(ax, 'Sum of squared errors');
title(ax, 'LS cost convergence'); grid(ax, 'on');
ax = nexttile(tl, 6); axis(ax, 'off');
trueVar = var(d.r.xb - d.xc);
txt = {'\bfFitted parameters\rm', '', ...
       sprintf('%-4s %8s %8s %8s', 'wave', '\alpha (mV)', 'b (rad)', '\theta (rad)')};
for i = 1:n
    txt{end+1} = sprintf('  %s    %7.3f   %6.3f   %6.3f', wn{i}, p(i), p(n+i), p(2*n+i)); %#ok<SAGROW>
end
txt = [txt, {'', '\bfML noise variance\rm (Kalman R)', ...
       sprintf('\\sigma^2_{MLE} = %.4f mV^2', d.r.sigma2), ...
       sprintf('true value = %.4f mV^2', trueVar)}];
text(ax, 0, 1, txt, 'VerticalAlignment', 'top', 'FontSize', 15, 'FontName', 'Consolas');
grab_frame(f, fullfile(figDir, 'fig03_ls_gaussian_fit.png')); close(f);

%% ===== Fig 4: denoising comparison at one noise level =====
[~, f] = fig_style(true);
tl = tiledlayout(f, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
k = tw(20, 24);
ax = nexttile(tl, 1, [1 2]);
plot(ax, t(k), d.y(k), 'Color', C.noisy, 'LineWidth', 1.2); hold(ax, 'on');
plot(ax, t(k), d.xc(k), 'k', 'LineWidth', 2);
title(ax, sprintf('Noisy input (grey, MA + BW at %d dB)   vs   clean reference ECG (black)', d.snr));
ylabel(ax, 'mV'); grid(ax, 'on');
est = {d.bp, d.wt, d.r.xEKF, d.r.xEKS};
for j = 1:4
    ax = nexttile(tl);
    plot(ax, t(k), d.xc(k), 'k', 'LineWidth', 1.6); hold(ax, 'on');
    plot(ax, t(k), est{j}(k), 'Color', C.mcol{j}, 'LineWidth', 2.2);
    title(ax, sprintf('%s:   SNR_{out} %.1f dB,  RMSE %.3f mV', C.short{j}, ...
          d.m(j).snrOut, d.m(j).rmse), 'Color', C.mcol{j}*0.85, 'FontSize', 15);
    ylim(ax, [-0.8 1.6]); grid(ax, 'on'); ylabel(ax, 'mV');
    if j > 2, xlabel(ax, 'Time (s)'); end
end
grab_frame(f, fullfile(figDir, 'fig04_method_comparison.png')); close(f);

%% ===== Fig 5: inside the Kalman filter - uncertainty and gain =====
[~, f] = fig_style(true);
tl = tiledlayout(f, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
k = tw(20, 23);  tk = t(k);
ax = nexttile(tl); hold(ax, 'on');
sd = 2*sqrt(d.r.ekf.Pzf(k));
fill(ax, [tk; flipud(tk)], [d.r.xEKF(k) + sd; flipud(d.r.xEKF(k) - sd)], C.ekf, ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(ax, tk, d.r.xb(k), 'Color', [C.noisy 0.6], 'LineWidth', 1);
plot(ax, tk, d.xc(k), 'k', 'LineWidth', 1.6);
plot(ax, tk, d.r.xEKF(k), 'Color', C.ekf, 'LineWidth', 2.4);
legend(ax, {'EKF \pm2\sigma (posterior covariance)', 'Noisy input', 'Clean', 'EKF estimate'}, ...
       'Location', 'eastoutside');
title(ax, 'EKF estimate with its own uncertainty'); ylabel(ax, 'mV'); grid(ax, 'on');
ax = nexttile(tl);
plot(ax, tk, d.r.ekf.Kz(k), 'Color', C.eks, 'LineWidth', 2.4);
title(ax, 'Kalman gain K_z: how much the filter trusts each new noisy sample');
xlabel(ax, 'Time (s)'); ylabel(ax, 'K_{z}'); grid(ax, 'on');
legend(ax, {'K_z(t)'}, 'Location', 'eastoutside');
grab_frame(f, fullfile(figDir, 'fig05_kalman_gain_uncertainty.png')); close(f);

%% ===== Fig 6: performance vs input SNR (main quantitative result) =====
[~, f] = fig_style(true);
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
tn = {'Muscle artifact (MA)', 'Baseline wander (BW)', 'MA + BW'};
for a = 1:3
    ax = nexttile(tl); hold(ax, 'on');
    bl = R.S.snrBL(:,a,:,:);  bl = reshape(permute(bl, [3 1 2 4]), numel(R.snrs), []);
    hb = plot(ax, R.snrs, mean(bl, 2), ':', 'Color', [0.35 0.35 0.35], 'LineWidth', 2.5);
    h = gobjects(1,4);
    for j = 1:4
        v = R.S.snrOut(:,a,:,:,j);  v = reshape(permute(v, [3 1 2 4]), numel(R.snrs), []);
        mu = mean(v, 2);  sdv = std(v, 0, 2);
        fill(ax, [R.snrs fliplr(R.snrs)], [mu - sdv; flipud(mu + sdv)]', C.mcol{j}, ...
             'FaceAlpha', 0.15, 'EdgeColor', 'none');
        h(j) = plot(ax, R.snrs, mu, 'o-', 'Color', C.mcol{j}, 'MarkerFaceColor', C.mcol{j}, ...
                    'LineWidth', 2.5, 'MarkerSize', 7);
    end
    title(ax, tn{a}); xlabel(ax, 'Input SNR (dB, nst)'); grid(ax, 'on');
    xlim(ax, [-5 20]); xticks(ax, -5:5:20);
    if a == 1, ylabel(ax, 'Output SNR (dB)'); end
end
lg = legend(ax, [h hb], [C.names, {'No denoising (baseline removed only)'}], 'Orientation', 'horizontal');
lg.Layout.Tile = 'south'; lg.FontSize = 15;
title(tl, sprintf('Output SNR vs noise level  (mean \\pm 1 SD, %d runs per point: records 100 & 103 \\times %d noise segments)', ...
      numel(R.rec)*R.nTrial, R.nTrial), 'FontSize', 18, 'FontWeight', 'bold');
grab_frame(f, fullfile(figDir, 'fig06_snr_sweep.png')); close(f);

%% ===== Fig 7: limitation - EKF gain depends on R-peak detection =====
[~, f] = fig_style(true);
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nc = {C.wt, C.ekf, C.eks};
ax = nexttile(tl); hold(ax, 'on');
for a = 1:3
    se = squeeze(mean(R.S.sens(:,a,:,:), [1 4]));  pp = squeeze(mean(R.S.ppv(:,a,:,:), [1 4]));
    plot(ax, R.snrs, 100*se, 'o-', 'Color', nc{a}, 'MarkerFaceColor', nc{a}, 'LineWidth', 2.5);
    plot(ax, R.snrs, 100*pp, 's--', 'Color', nc{a}, 'LineWidth', 2);
end
legend(ax, {'MA  sensitivity', 'MA  precision', 'BW  sensitivity', 'BW  precision', ...
            'MA+BW  sensitivity', 'MA+BW  precision'}, 'Location', 'southeast');
xlabel(ax, 'Input SNR (dB, nst)'); ylabel(ax, '%'); ylim(ax, [0 105]); grid(ax, 'on');
title(ax, 'R-peak detection on the noisy signal');
ax = nexttile(tl); hold(ax, 'on');
for a = 1:3
    g = R.S.snrOut(:,a,:,:,4) - R.S.snrOut(:,a,:,:,1);
    plot(ax, R.snrs, squeeze(mean(g, [1 4])), 'o-', 'Color', nc{a}, 'MarkerFaceColor', nc{a}, 'LineWidth', 2.5);
end
yline(ax, 0, 'k-', 'LineWidth', 1.2);
legend(ax, {'MA', 'BW', 'MA+BW'}, 'Location', 'northeast');
xlabel(ax, 'Input SNR (dB, nst)'); ylabel(ax, 'dB'); grid(ax, 'on');
title(ax, 'EKS gain over band-pass (\Delta output SNR)');
grab_frame(f, fullfile(figDir, 'fig07_limitation_rpeaks.png')); close(f);

%% ===== Fig 8: metrics table at 5 dB =====
[~, f] = fig_style(true);
ax = axes(f, 'Position', [0 0 1 1]); axis(ax, 'off'); hold(ax, 'on');
xlim(ax, [0 1]); ylim(ax, [0 1]);
b = find(R.snrs == 5);
text(ax, 0.5, 0.93, 'Quantitative comparison at input SNR = 5 dB', ...
     'FontSize', 26, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(ax, 0.5, 0.875, 'mean of 10 test runs (records 100 & 103, 5 noise segments each)', ...
     'FontSize', 16, 'Color', [0.35 0.35 0.35], 'HorizontalAlignment', 'center');
x0 = 0.04;  x1 = 0.96;                     % table extent
grpX = [0.44 0.66 0.86];                   % centre of each noise-type group
subX = [-0.07 0 0.065];                    % SNR_out, RMSE, r inside a group
grpName = {'Muscle artifact', 'Baseline wander', 'MA + BW'};
sub = {{'SNR_{out}', '(dB)'}, {'RMSE', '(mV)'}, {'r', ''}};
for a = 1:3
    text(ax, grpX(a), 0.78, grpName{a}, 'FontSize', 20, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    plot(ax, grpX(a) + [-0.1 0.1], [0.75 0.75], 'k-', 'LineWidth', 1.2);
    for c = 1:3
        text(ax, grpX(a) + subX(c), 0.722, sub{c}, 'FontSize', 13, 'Color', [0.3 0.3 0.3], ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
text(ax, x0 + 0.01, 0.715, 'Method', 'FontSize', 18, 'FontWeight', 'bold');
plot(ax, [x0 x1], [0.68 0.68], 'k-', 'LineWidth', 1.5);
subName = {'Butterworth 0.5-40 Hz', 'sym4, empirical Bayes', 'forward filter (causal)', 'RTS smoother (offline)'};
for j = 1:4
    y0 = 0.60 - (j-1)*0.135;
    rectangle(ax, 'Position', [x0 y0-0.06 x1-x0 0.12], 'FaceColor', [C.mcol{j} 0.12], 'EdgeColor', 'none');
    text(ax, x0 + 0.01, y0 + 0.018, C.short{j}, 'FontSize', 21, 'FontWeight', 'bold', 'Color', C.mcol{j}*0.85);
    text(ax, x0 + 0.01, y0 - 0.03, subName{j}, 'FontSize', 13, 'Color', [0.35 0.35 0.35]);
    for a = 1:3
        vals = {sprintf('%.1f', mean(R.S.snrOut(:,a,b,:,j), 'all')), ...
                sprintf('%.3f', mean(R.S.rmse(:,a,b,:,j), 'all')), ...
                sprintf('%.2f', mean(R.S.corr(:,a,b,:,j), 'all'))};
        for c = 1:3
            text(ax, grpX(a) + subX(c), y0, vals{c}, 'FontSize', 20, 'HorizontalAlignment', 'center', ...
                 'FontWeight', ifelse(j == 4, 'bold', 'normal'));
        end
    end
end
text(ax, 0.5, 0.06, ['SNR_{out}: higher is better      RMSE: lower is better      ' ...
     'r: correlation with the clean ECG (1 = perfect)'], 'FontSize', 15, ...
     'Color', [0.3 0.3 0.3], 'HorizontalAlignment', 'center');
grab_frame(f, fullfile(figDir, 'fig08_metrics_table.png')); close(f);

%% ===== Fig 9: real noisy recording (MIT-BIH 105, 20:00-21:00) - no ground truth =====
[~, f] = fig_style(true);
rl = R.real(1);  tr = rl.t0 + (0:R.L-1)'/fs;
w0 = rl.t0 + 20;  k = tr >= w0 & tr < w0 + 6;
beatCodes = [1:13 25 30 31 34 35 37 38];
ann = rl.annSamp(ismember(rl.annType, beatCodes));  ann = ann(k(ann));
det = rl.r.pk(k(rl.r.pk));
tl = tiledlayout(f, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); hold(ax, 'on');
plot(ax, tr(k), rl.y(k), 'Color', C.noisy*0.7, 'LineWidth', 1.4);
yl = ylim(ax);  plot(ax, tr(ann), (yl(2) - 0.1)*ones(size(ann)), 'kv', 'MarkerSize', 11, 'MarkerFaceColor', 'k');
title(ax, 'Record 105 (20:20): real recording with strong natural noise   (\nabla = annotated beats)');
ylabel(ax, 'mV'); grid(ax, 'on');
ax = nexttile(tl); plot(ax, tr(k), rl.r.xb(k), 'Color', [C.noisy 0.6], 'LineWidth', 1); hold(ax, 'on');
plot(ax, tr(k), rl.bp(k), 'Color', C.bp, 'LineWidth', 2.2);
title(ax, 'Band-pass 0.5-40 Hz: keeps every beat, but also most of the noise'); ylabel(ax, 'mV'); grid(ax, 'on');
ax = nexttile(tl); plot(ax, tr(k), rl.r.xb(k), 'Color', [C.noisy 0.6], 'LineWidth', 1); hold(ax, 'on');
plot(ax, tr(k), rl.r.xEKS(k), 'Color', C.eks, 'LineWidth', 2.4);
yl = ylim(ax);  plot(ax, tr(det), (yl(2) - 0.15)*ones(size(det)), 'v', 'MarkerSize', 11, ...
                     'MarkerFaceColor', C.eks, 'MarkerEdgeColor', C.eks);
title(ax, 'EKS: clean beat shape, but beats are lost where R-peak detection fails   (\nabla = detected peaks)');
ylabel(ax, 'mV'); xlabel(ax, 'Time (s)'); grid(ax, 'on');
grab_frame(f, fullfile(figDir, 'fig09_real_record105_noisy.png')); close(f);

%% ===== Fig 10: limitation - abnormal beat (PVC) is suppressed by the model =====
[~, f] = fig_style(true);
rl = R.real(2);  tr = rl.t0 + (0:R.L-1)'/fs;
pvc = rl.annSamp(rl.annType == 5);                       % PVC annotations (code 5 = V)
pvc = pvc(pvc > 4*fs & pvc < R.L - 4*fs);                % well inside the segment
w = round(0.15*fs);                                      % pick the PVC the EKS distorts most
err = arrayfun(@(p) sqrt(mean((rl.r.xb(p-w:p+w) - rl.r.xEKS(p-w:p+w)).^2)), pvc);
[~, iw] = max(err);
tc  = tr(pvc(iw));  k = tr >= tc - 2.2 & tr < tc + 2.2;
tl = tiledlayout(f, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
yr = [min(rl.r.xb(k)) - 0.3, max(rl.r.xb(k)) + 0.6];
lbl = {'Band-pass keeps the PVC', 'EKS pulls the PVC towards the normal-beat template (shape lost)'};
est = {rl.bp, rl.r.xEKS};  col = {C.bp, C.eks};  nm = {'Band-pass', 'EKS'};
for i = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    patch(ax, tc + [-0.12 0.12 0.12 -0.12], yr([1 1 2 2]), [1 0.85 0.85], 'EdgeColor', 'none', ...
          'HandleVisibility', 'off');
    plot(ax, tr(k), rl.r.xb(k), 'Color', C.noisy*0.6, 'LineWidth', 2);
    plot(ax, tr(k), est{i}(k), 'Color', col{i}, 'LineWidth', 2.6);
    text(ax, tc, yr(2) - 0.08, 'PVC', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.7 0 0], ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    ylim(ax, yr); xlim(ax, [tc - 2.2, tc + 2.2]); grid(ax, 'on'); ylabel(ax, 'mV');
    legend(ax, {'Recording (baseline removed)', nm{i}}, 'Location', 'northeast');
    title(ax, lbl{i});
end
xlabel(ax, 'Time (s)');
title(tl, 'Limitation: the model assumes every beat has the same shape (MIT-BIH record 105)', ...
      'FontSize', 20, 'FontWeight', 'bold');
grab_frame(f, fullfile(figDir, 'fig10_limitation_pvc.png')); close(f);

%% ===== Fig 11: generalisation - the unseen patient (record 103) =====
[~, f] = fig_style(true);
q2 = find([R.rec.id] == 103);
y103 = R.rec(q2).xc + make_noise(ds, 'ma+bw', 5, R.L, R.offsets(1), R.rec(q2).P);
r103 = ekf_pipeline(y103, fs, R.opt);
m103 = eval_metrics(R.rec(q2).xc, y103, r103.xEKS, fs);
bp103 = denoise_bandpass(r103.xb, fs);  mb103 = eval_metrics(R.rec(q2).xc, y103, bp103, fs);
tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl, 1, [1 3]); hold(ax, 'on');
k = tw(20, 24);
plot(ax, t(k), r103.xb(k), 'Color', C.noisy, 'LineWidth', 1.1);
plot(ax, t(k), R.rec(q2).xc(k), 'k', 'LineWidth', 1.8);
plot(ax, t(k), r103.xEKS(k), 'Color', C.eks, 'LineWidth', 2.4);
title(ax, sprintf('Record 103, MA + BW at 5 dB:  EKS %.1f dB vs band-pass %.1f dB   (grey noisy, black clean, orange EKS)', ...
      m103.snrOut, mb103.snrOut), 'FontSize', 15);
ylabel(ax, 'mV'); xlabel(ax, 'Time (s)'); grid(ax, 'on');
tn = {'Muscle artifact (MA)', 'Baseline wander (BW)', 'MA + BW'};
for a = 1:3
    ax = nexttile(tl); hold(ax, 'on');
    h = gobjects(1,4);  n = 0;
    for j = [1 4]                                   % band-pass and EKS
        for q = 1:numel(R.rec)
            v = squeeze(mean(R.S.snrOut(q,a,:,:,j), 4));
            ls = '-';  mk = 'o';
            if R.rec(q).id == 103, ls = '--'; mk = 's'; end
            n = n + 1;
            h(n) = plot(ax, R.snrs, v, [mk ls], 'Color', C.mcol{j}, 'MarkerFaceColor', C.mcol{j}, ...
                        'LineWidth', 2.4, 'MarkerSize', 6);
        end
    end
    title(ax, tn{a}); xlabel(ax, 'Input SNR (dB, nst)'); grid(ax, 'on'); xlim(ax, [-5 20]);
    if a == 1, ylabel(ax, 'Output SNR (dB)'); end
end
lg = legend(ax, h, {'Band-pass, record 100', 'Band-pass, record 103 (unseen)', ...
                    'EKS, record 100', 'EKS, record 103 (unseen)'}, 'Orientation', 'horizontal');
lg.Layout.Tile = 'south'; lg.FontSize = 15;
title(tl, 'Same performance on a patient never used for tuning (record 103)', ...
      'FontSize', 18, 'FontWeight', 'bold');
grab_frame(f, fullfile(figDir, 'fig11_unseen_patient_103.png')); close(f);

%% ===== Fig 12: how the ML noise variance is estimated =====
[~, f] = fig_style(true);
res = d.r.xb - gauss_ecg_model(d.r.p, d.r.phi);           % residual = data - LS beat model
sig = sqrt(d.r.sigma2);  trueSig = std(d.r.xb - d.xc);
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');
ax = nexttile(tl); hold(ax, 'on');
k = tw(20, 23);
plot(ax, t(k), d.r.xb(k), 'Color', C.noisy, 'LineWidth', 1.2);
plot(ax, t(k), gauss_ecg_model(d.r.p, d.r.phi(k)), 'Color', C.eks, 'LineWidth', 2.4);
plot(ax, t(k), res(k) - 1.6, 'Color', C.ekf, 'LineWidth', 1.4);
text(ax, 20.05, -0.82, 'residual e = data - model', 'FontSize', 16, 'Color', C.ekf, 'FontWeight', 'bold');
legend(ax, {'Noisy ECG (data)', 'LS beat model z(\theta)', 'Residual (shifted down)'}, 'Location', 'northeast');
title(ax, '1) Subtract the LS beat model from the data');
xlabel(ax, 'Time (s)'); ylabel(ax, 'mV'); grid(ax, 'on'); ylim(ax, [-2.4 1.8]);
ax = nexttile(tl); hold(ax, 'on');
edges = linspace(-1, 1, 61);
histogram(ax, res, edges, 'Normalization', 'pdf', 'FaceColor', C.ekf, 'FaceAlpha', 0.35, 'EdgeColor', 'none');
xx = linspace(-1, 1, 400);
plot(ax, xx, exp(-xx.^2/(2*sig^2))/(sqrt(2*pi)*sig), 'Color', C.eks, 'LineWidth', 3);
title(ax, '2) Fit a Gaussian N(0, \sigma^2) by maximum likelihood');
xlabel(ax, 'Residual (mV)'); ylabel(ax, 'Probability density'); grid(ax, 'on'); xlim(ax, [-1 1]);
legend(ax, {'Residual histogram', sprintf('ML Gaussian, \\sigma^2 = %.4f mV^2', d.r.sigma2)}, 'Location', 'northeast');
text(ax, -0.97, 0.82*max(ylim(ax)), {'\sigma^2_{ML} = mean(e^2)', sprintf('= %.4f mV^2', d.r.sigma2), ...
     sprintf('true: %.4f mV^2', trueSig^2)}, 'FontSize', 16, 'VerticalAlignment', 'top');
grab_frame(f, fullfile(figDir, 'fig12_mle_noise_variance.png')); close(f);

fprintf('Figures saved to %s\n', figDir);

function out = ifelse(c, a, b)
    if c, out = a; else, out = b; end
end
