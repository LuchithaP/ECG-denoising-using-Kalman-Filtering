function [S, N] = nst_power(ecg, annSamp, annType, fs, noise)
%NST_POWER  Signal and noise power exactly as defined by PhysioNet's nst.
%   S : QRS "power" = (trimmed-mean peak-to-peak QRS amplitude)^2 / 8,
%       measured over the first 300 normal beats, window +-50 ms.
%   N : noise power = (trimmed-mean RMS of 1-s chunks, about the chunk mean)^2,
%       measured over the first 300 s of the noise record.
%   In both cases the largest 5 % and smallest 5 % of measurements are dropped.
%   Reference: https://physionet.org/physiotools/wag/nst-1.htm

    S = []; N = [];
    if ~isempty(ecg)
        w   = round(0.05*fs);
        idx = annSamp(annType == 1);                 % normal beats only
        idx = idx(idx > w & idx <= numel(ecg) - w);
        idx = idx(1:min(300, numel(idx)));
        pp  = arrayfun(@(k) max(ecg(k-w:k+w)) - min(ecg(k-w:k+w)), idx);
        S   = trimmed_mean(pp)^2 / 8;
    end
    if nargin > 4 && ~isempty(noise)
        nCh   = min(300, floor(numel(noise)/fs));
        chunk = reshape(noise(1:nCh*fs), fs, nCh);
        n     = sqrt(mean((chunk - mean(chunk)).^2));  % RMS about chunk mean
        N     = trimmed_mean(n)^2;
    end
end

function m = trimmed_mean(v)
    v = sort(v(:));
    k = floor(0.05*numel(v));
    m = mean(v(k+1:end-k));
end
