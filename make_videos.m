%% MAKE_VIDEOS  Short MP4 clips (1920x1080, 30 fps) for the presentation video.
%  vid01_ekf_live.mp4  : the EKF tracking a noisy ECG in real time, with the
%                        state (phase, amplitude) shown on a "phase wheel"
%  vid02_snr_sweep.mp4 : noise level increasing from 20 dB to -5 dB; the
%                        filtered output and the output-SNR curves update live
%  Requires results/experiment_results.mat.  
%
%  Robust rendering: every frame is first saved as a PNG in
%  results/videos/frames_*/ and frames that already exist are skipped, so if
%  MATLAB crashes during a long render simply run this script again - it
%  resumes where it stopped.  The MP4 files are assembled at the end.
%  Delete the frames_* folders to force a full re-render.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'));
ds = load(fullfile(here, 'data', 'ecg_dataset.mat'));
R  = load(fullfile(here, 'results', 'experiment_results.mat'));
vidDir = fullfile(here, 'results', 'videos');
C = fig_style();  fs = R.fs;  d = R.demo;  fps = 30;
fr1 = fullfile(vidDir, 'frames_vid01');  fr2 = fullfile(vidDir, 'frames_vid02');
if ~exist(fr1, 'dir'), mkdir(fr1); end
if ~exist(fr2, 'dir'), mkdir(fr2); end
frameFile = @(dir, n) fullfile(dir, sprintf('f%04d.png', n));

%% ================= Video 1: EKF in action =================
tStart = 20;  dur = 10;  win = 3;              % show 10 s of signal, 3 s window
step = round(fs/fps);                           % real-time speed
t  = (0:R.L-1)'/fs;
i0 = round(tStart*fs);  iEnd = round((tStart + dur)*fs);
idxList = i0:step:iEnd;
if ~all(arrayfun(@(n) exist(frameFile(fr1, n), 'file') == 2, 1:numel(idxList)))
    [~, f] = fig_style(true);
    axT = axes(f, 'Position', [0.07 0.19 0.56 0.62]); hold(axT, 'on'); grid(axT, 'on');
    hN = plot(axT, nan, nan, 'Color', C.noisy, 'LineWidth', 1.2);
    hC = plot(axT, nan, nan, 'k', 'LineWidth', 1.6);
    hE = plot(axT, nan, nan, 'Color', C.ekf, 'LineWidth', 3);
    hD = plot(axT, nan, nan, 'o', 'MarkerSize', 12, 'MarkerFaceColor', C.ekf, 'MarkerEdgeColor', 'w');
    ylim(axT, [-1.2 2.0]); xlabel(axT, 'Time (s)'); ylabel(axT, 'Amplitude (mV)');
    legend(axT, {'Noisy input (MA+BW, 5 dB)', 'Clean ECG', 'EKF estimate'}, ...
           'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', 15);
    axP = axes(f, 'Position', [0.67 0.14 0.31 0.72]); hold(axP, 'on'); axis(axP, 'equal', 'off');
    th = linspace(-pi, pi, 500)';  zT = gauss_ecg_model(d.r.p, th);
    rad = @(z) 1 + 0.45*z;                         % amplitude -> radius
    plot(axP, cos(th), sin(th), ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    plot(axP, rad(zT).*cos(th), rad(zT).*sin(th), 'Color', [C.eks 0.6], 'LineWidth', 3);
    hTr = plot(axP, nan, nan, 'Color', C.ekf, 'LineWidth', 2.5);
    hPt = plot(axP, nan, nan, 'o', 'MarkerSize', 16, 'MarkerFaceColor', C.ekf, 'MarkerEdgeColor', 'w');
    hAr = plot(axP, nan, nan, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
    text(axP, 0, 1.95, 'State space: phase wheel', 'FontSize', 20, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text(axP, 0, -1.95, {'angle = cardiac phase \theta', 'radius = ECG amplitude z'}, 'FontSize', 15, 'HorizontalAlignment', 'center');
    text(axP, 1.05, 0, '  R', 'FontSize', 18, 'FontWeight', 'bold', 'Color', C.eks);
    xlim(axP, [-2 2]); ylim(axP, [-2.1 2.1]);
    hTxt = annotation(f, 'textbox', [0.05 0.005 0.9 0.06], 'String', '', 'FontSize', 18, ...
                      'EdgeColor', 'none', 'HorizontalAlignment', 'left');
    annotation(f, 'textbox', [0.05 0.92 0.9 0.07], 'String', ...
        'Extended Kalman Filter tracking a noisy ECG in real time', 'FontSize', 26, ...
        'FontWeight', 'bold', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    for n = 1:numel(idxList)
        if exist(frameFile(fr1, n), 'file'), continue; end
        i = idxList(n);
        a = max(1, i - round(win*fs));  k = a:i;
        set(hN, 'XData', t(k), 'YData', d.r.xb(k));                    % baseline-removed input
        set(hC, 'XData', t(k), 'YData', d.xc(k));
        set(hE, 'XData', t(k), 'YData', d.r.xEKF(k));
        set(hD, 'XData', t(i), 'YData', d.r.xEKF(i));
        xlim(axT, [t(i) - win, t(i) + 0.15]);
        kk = max(1, i - round(0.8*fs)):i;           % recent trajectory on the wheel
        th_ = d.r.ekf.thf(kk);  z_ = d.r.xEKF(kk);
        set(hTr, 'XData', rad(z_).*cos(th_), 'YData', rad(z_).*sin(th_));
        set(hPt, 'XData', rad(z_(end))*cos(th_(end)), 'YData', rad(z_(end))*sin(th_(end)));
        set(hAr, 'XData', [0 rad(z_(end))*cos(th_(end))], 'YData', [0 rad(z_(end))*sin(th_(end))]);
        set(hTxt, 'String', sprintf(['t = %5.2f s     heart rate = %3.0f bpm     ' ...
            'Kalman gain K_z = %.3f     posterior SD = %.3f mV'], t(i), 60*d.r.omega(i)/(2*pi), ...
            d.r.ekf.Kz(i), sqrt(d.r.ekf.Pzf(i))));
        grab_frame(f, frameFile(fr1, n));
    end
    close(f);
end
assemble_video(fullfile(vidDir, 'vid01_ekf_live.mp4'), fr1, numel(idxList), fps, 1, fps);
fprintf('Saved vid01_ekf_live.mp4\n');

%% ================= Video 2: noise-level sweep =================
snrV = 20:-0.5:-5;  holdF = 6;                  % 0.2 s per noise level
q = 1;  xc = R.rec(q).xc;  P = R.rec(q).P;  off = R.offsets(1);
Ls = 30*fs;  xs = xc(1:Ls);  ts = (0:Ls-1)'/fs;  kw = ts >= 10 & ts < 14;
cacheFile = fullfile(fr2, 'sweep_cache.mat');
if exist(cacheFile, 'file')
    load(cacheFile, 'res', 'segN', 'segB', 'segS');
else                                            % compute all levels first (fast)
    res = nan(numel(snrV), 4);  segN = nan(nnz(kw), numel(snrV));  segB = segN;  segS = segN;
    for s = 1:numel(snrV)
        y = xs + make_noise(ds, 'ma+bw', snrV(s), Ls, off, P);
        r = ekf_pipeline(y, fs, R.opt);
        est = {denoise_bandpass(r.xb, fs), denoise_wavelet(r.xb), r.xEKF, r.xEKS};
        for j = 1:4, m = eval_metrics(xs, y, est{j}, fs); res(s,j) = m.snrOut; end
        segN(:,s) = r.xb(kw);  segB(:,s) = est{1}(kw);  segS(:,s) = est{4}(kw);
    end
    save(cacheFile, 'res', 'segN', 'segB', 'segS', 'snrV');
end
if ~all(arrayfun(@(n) exist(frameFile(fr2, n), 'file') == 2, 1:numel(snrV)))
    [~, f] = fig_style(true);
    tl = tiledlayout(f, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    axS = nexttile(tl, 1, [1 3]); hold(axS, 'on'); grid(axS, 'on');
    hN = plot(axS, ts(kw), nan(nnz(kw),1), 'Color', C.noisy, 'LineWidth', 1.1);
    plot(axS, ts(kw), xs(kw), 'k', 'LineWidth', 1.6);
    hB = plot(axS, ts(kw), nan(nnz(kw),1), 'Color', C.bp, 'LineWidth', 2);
    hS = plot(axS, ts(kw), nan(nnz(kw),1), 'Color', C.eks, 'LineWidth', 2.6);
    ylim(axS, [-1.5 2.2]); ylabel(axS, 'mV'); xlabel(axS, 'Time (s)');
    legend(axS, {'Noisy input', 'Clean', 'Band-pass', 'EKS'}, 'Location', 'eastoutside', 'FontSize', 15);
    hT = title(axS, '', 'FontSize', 22);
    axR = nexttile(tl, 4, [1 3]); hold(axR, 'on'); grid(axR, 'on');
    hl = gobjects(1,4);
    for j = 1:4
        hl(j) = plot(axR, nan, nan, 'o-', 'Color', C.mcol{j}, 'MarkerFaceColor', C.mcol{j}, ...
                     'LineWidth', 2.5, 'MarkerSize', 5);
    end
    set(axR, 'XDir', 'reverse'); xlim(axR, [-5 20]); ylim(axR, [-8 26]);
    xlabel(axR, 'Input SNR (dB, nst)   \leftarrow  noise increases'); ylabel(axR, 'Output SNR (dB)');
    legend(axR, C.names, 'Location', 'eastoutside', 'FontSize', 15);
    for s = 1:numel(snrV)
        if exist(frameFile(fr2, s), 'file'), continue; end
        set(hN, 'YData', segN(:,s));  set(hB, 'YData', segB(:,s));  set(hS, 'YData', segS(:,s));
        set(hT, 'String', sprintf('Input SNR %5.1f dB   |   output: EKS %5.1f dB,  band-pass %5.1f dB', ...
            snrV(s), res(s,4), res(s,1)));
        for j = 1:4, set(hl(j), 'XData', snrV(1:s), 'YData', res(1:s,j)); end
        grab_frame(f, frameFile(fr2, s));
    end
    close(f);
end
assemble_video(fullfile(vidDir, 'vid02_snr_sweep.mp4'), fr2, numel(snrV), fps, holdF, 2*fps);
fprintf('Saved vid02_snr_sweep.mp4\n');

function assemble_video(outFile, frameDir, nFrames, fps, repeatEach, holdLast)
%ASSEMBLE_VIDEO  Build an MP4 from numbered PNG frames.
%   repeatEach : how many video frames each PNG is shown for
%   holdLast   : extra video frames to hold the last PNG
    vw = VideoWriter(outFile, 'MPEG-4');
    vw.FrameRate = fps;  vw.Quality = 95;  open(vw);
    for n = 1:nFrames
        im = imread(fullfile(frameDir, sprintf('f%04d.png', n)));
        for r = 1:repeatEach, writeVideo(vw, im); end
    end
    for r = 1:holdLast, writeVideo(vw, im); end
    close(vw);
end
