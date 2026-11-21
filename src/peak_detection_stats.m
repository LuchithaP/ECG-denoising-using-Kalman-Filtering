function [sens, ppv] = peak_detection_stats(pkDet, pkRef, fs)
%PEAK_DETECTION_STATS  R-peak detection sensitivity and positive predictivity.
%   A detection is a true positive if it lies within 50 ms of a reference
%   beat annotation (each reference matched at most once).
    tol = round(0.05*fs);
    used = false(size(pkRef));  tp = 0;
    for i = 1:numel(pkDet)
        [dmin, j] = min(abs(pkRef - pkDet(i)) + 1e9*used);
        if dmin <= tol, tp = tp + 1; used(j) = true; end
    end
    sens = tp / max(numel(pkRef), 1);
    ppv  = tp / max(numel(pkDet), 1);
end
