function sp = data_split(ds)
%DATA_SPLIT  Single place that defines the TRAIN (tuning) and TEST data.
%
%  TRAIN (used ONLY to tune the EKF process-noise hyper-parameters)
%    - clean ECG : record 100, t = 1:00-2:00 (60 s)
%    - noise     : one fixed NSTDB segment starting at sample 400000
%
%  TEST (never seen during tuning)
%    - clean ECG : record 100, t = 5:00-6:00   and   record 103, t = 5:00-6:00
%    - noise     : 5 random NSTDB segments (seed 2026) that do not overlap
%                  the training noise segment
%
%  NOTE: the LS wave parameters and the ML noise variance are NOT trained;
%  they are re-estimated from every noisy test recording itself (no clean
%  signal or labels are used), so nothing leaks from the test set.

    fs = ds.fs;  L = 60*fs;
    sp.fs = fs;  sp.L = L;

    sp.train.rec    = 100;
    sp.train.start  = 60*fs;          % sample offset (1:00)
    sp.train.noiseOffset = 400000;

    sp.test.rec   = [100 103];
    sp.test.start = [300 300]*fs;     % 5:00 in both records
    nTrial = 5;
    rs = RandStream('mt19937ar', 'Seed', 2026);
    nNoise = numel(ds.noise_ma);
    off = zeros(nTrial,1);
    for t = 1:nTrial
        o = sp.train.noiseOffset;
        while abs(o - sp.train.noiseOffset) < 2*L, o = randi(rs, nNoise - L); end
        off(t) = o;
    end
    sp.test.noiseOffsets = off;
end
