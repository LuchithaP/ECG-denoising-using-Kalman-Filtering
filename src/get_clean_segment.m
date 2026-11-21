function [xc, pkRef, S] = get_clean_segment(ds, rec, start, L)
%GET_CLEAN_SEGMENT  Ground-truth ECG segment, its reference beats and nst power.
%   xc    : baseline-removed clean ECG (L samples) = ground truth
%   pkRef : annotated beat locations inside the segment (for R-peak scoring)
%   S     : nst signal power of the record (first 300 normal beats)
    ecg  = ds.(sprintf('ecg%d', rec));
    aS   = ds.(sprintf('ann%d_samp', rec));
    aT   = ds.(sprintf('ann%d_type', rec));
    xc   = remove_baseline(ecg(start+1:start+L), ds.fs);
    beatCodes = [1:13 25 30 31 34 35 37 38];      % WFDB beat annotation codes
    k = ismember(aT, beatCodes) & aS > start & aS <= start + L;
    pkRef = aS(k) - start;
    S = nst_power(ecg, aS, aT, ds.fs);
end
