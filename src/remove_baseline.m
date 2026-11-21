function [y, base] = remove_baseline(x, fs)
%REMOVE_BASELINE  Two-stage median filter baseline estimation (200 ms, 600 ms).
%   The 200 ms median removes QRS complexes and P waves, the 600 ms median
%   removes T waves; what remains is the baseline wander, which is subtracted.
%   This common front-end is applied before every denoising method.
    w1 = 2*floor(0.2*fs/2) + 1;
    w2 = 2*floor(0.6*fs/2) + 1;
    base = medfilt1(medfilt1(x(:), w1, 'truncate'), w2, 'truncate');
    y = x(:) - base;
end
