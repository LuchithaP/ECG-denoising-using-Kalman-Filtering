function [C, f] = fig_style(makeFig)
%FIG_STYLE  Common colours and a 1920x1080 (16:9) figure for slides/video.
%   Colour-blind-safe palette (Okabe-Ito).  Large fonts and thick lines so
%   that every plot stays readable when shown full-screen in a video.
    C.clean = [0 0 0];
    C.noisy = [0.62 0.62 0.62];
    C.bp    = [0 158 115]/255;     % band-pass  (green)
    C.wt    = [204 121 167]/255;   % wavelet    (pink)
    C.ekf   = [0 114 178]/255;     % EKF        (blue)
    C.eks   = [213 94 0]/255;      % EKS        (orange)
    C.waves = [0.80 0.47 0.65; 0.34 0.71 0.91; 0.84 0.37 0.00; 0.00 0.62 0.45; 0.90 0.62 0.00];
    C.names = {'Band-pass 0.5-40 Hz', 'Wavelet (sym4, Bayes)', 'EKF', 'EKS (RTS smoother)'};
    C.short = {'Band-pass', 'Wavelet', 'EKF', 'EKS'};
    C.mcol  = {C.bp, C.wt, C.ekf, C.eks};
    f = [];
    if nargin > 0 && makeFig
        f = figure('Visible', 'off', 'Units', 'pixels', 'Position', [0 0 1920 1080], ...
                   'Color', 'w', 'PaperUnits', 'inches', ...
                   'PaperPosition', [0 0 12.8 7.2], 'PaperPositionMode', 'manual');
        set(f, 'DefaultAxesFontSize', 15, 'DefaultAxesFontName', 'Arial', ...
               'DefaultTextFontName', 'Arial', 'DefaultAxesLineWidth', 1.2, ...
               'DefaultLineLineWidth', 2, 'DefaultAxesTitleFontSizeMultiplier', 1.15, ...
               'DefaultAxesBox', 'on', 'DefaultLegendFontSize', 14);
    end
end
