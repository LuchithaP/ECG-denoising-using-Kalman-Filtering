function pk = detect_rpeaks(x, fs)
%DETECT_RPEAKS  Simple Pan-Tompkins style R-peak detector (works on noisy ECG).
% 1) band-pass 5-15 Hz (QRS energy band)  
% 2) derivative, squaring
% 3) 150 ms moving-window integration      
% 4) peak picking with a refractory period of 250 ms and a robust threshold
% 5) refine each detection to the largest |amplitude| in the 5-30 Hz signal.

    x = x(:);
    [b, a] = butter(3, [5 15]/(fs/2), 'bandpass');
    xb = filtfilt(b, a, x);
    e  = movmean(gradient(xb).^2, round(0.15*fs));   % integrated energy

    thr = 0.3 * prctile(e, 98);
    [~, loc] = findpeaks(e, 'MinPeakHeight', thr, 'MinPeakDistance', round(0.25*fs));

    % Refine to the R-wave extremum of a slightly wider band-passed signal
    [b2, a2] = butter(3, [5 30]/(fs/2), 'bandpass');
    xr = filtfilt(b2, a2, x);
    pol = sign(median(xr(loc)));  if pol == 0, pol = 1; end    % R-wave polarity
    w = round(0.08*fs);
    pk = zeros(size(loc));
    for i = 1:numel(loc)
        r = max(1, loc(i)-w) : min(numel(x), loc(i)+w);
        [~, j] = max(pol * xr(r));
        pk(i) = r(j);
    end
    pk = unique(pk);
    % remove detections closer than 250 ms (keep the larger one)
    keep = true(size(pk));
    for i = 2:numel(pk)
        if pk(i) - pk(i-1) < 0.25*fs
            if pol*xr(pk(i)) > pol*xr(pk(i-1)), keep(i-1) = false; else, keep(i) = false; end
        end
    end
    pk = pk(keep);
end
