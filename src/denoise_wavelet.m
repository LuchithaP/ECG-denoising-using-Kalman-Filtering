function y = denoise_wavelet(x)
%DENOISE_WAVELET  Conventional baseline: wavelet shrinkage denoising.
%   sym4 wavelet, 5 decomposition levels, empirical-Bayes thresholds with a
%   level-dependent noise estimate (suits coloured noise such as EMG).
    y = wdenoise(x(:), 5, 'Wavelet', 'sym4', 'DenoisingMethod', 'Bayes', ...
                 'NoiseEstimate', 'LevelDependent');
end
