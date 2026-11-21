function y = denoise_bandpass(x, fs)
%DENOISE_BANDPASS  Conventional baseline: zero-phase 4th-order Butterworth
%   band-pass filter, 0.5-40 Hz (standard monitoring-ECG bandwidth).
    [b, a] = butter(4, [0.5 40]/(fs/2), 'bandpass');
    y = filtfilt(b, a, x(:));
end
