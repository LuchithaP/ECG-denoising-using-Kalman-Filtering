function [sig, fs, info] = read_wfdb212(recPath)
%READ_WFDB212  Read a PhysioNet WFDB record stored in format 212.
%   [sig, fs, info] = read_wfdb212('data/raw/100') reads 100.hea / 100.dat
%   and returns the signals in physical units (mV), one column per channel.


    % ---------- parse header (.hea) ----------
    fid = fopen([recPath '.hea'], 'r');
    assert(fid > 0, 'Cannot open header %s.hea', recPath);
    first = strsplit(strtrim(fgetl(fid)));
    nSig  = str2double(first{2});
    fs    = str2double(first{3});
    nSamp = str2double(first{4});
    gain = zeros(1,nSig); base = zeros(1,nSig); desc = cell(1,nSig);
    for c = 1:nSig
        f = strsplit(strtrim(fgetl(fid)));
        assert(strcmp(f{2},'212'), 'Only format 212 is supported');
        g = f{3};                                  % e.g. "200", "200(0)/mV"
        gain(c) = str2double(regexp(g, '^[\d\.]+', 'match', 'once'));
        if gain(c) == 0, gain(c) = 200; end        % WFDB default gain
        adcZero = str2double(f{5});
        bTok = regexp(g, '\((-?\d+)\)', 'tokens', 'once');   % optional "(baseline)"
        if ~isempty(bTok), base(c) = str2double(bTok{1}); else, base(c) = adcZero; end
        desc{c} = strjoin(f(9:end), ' ');
    end
    fclose(fid);

    % ---------- read samples (.dat) ----------
    fid = fopen([recPath '.dat'], 'r');
    raw = fread(fid, inf, 'uint8=>uint16');
    fclose(fid);
    raw = raw(1:3*floor(numel(raw)/3));
    b = reshape(raw, 3, []);
    s1 = bitor(b(1,:), bitshift(bitand(b(2,:), 15), 8));      % low nibble of byte 2
    s2 = bitor(b(3,:), bitshift(bitand(b(2,:), 240), 4));     % high nibble of byte 2
    d  = double([s1; s2]);
    d(d > 2047) = d(d > 2047) - 4096;                         % 12-bit two's complement
    d  = reshape(d, nSig, []).';                              % samples x channels
    d  = d(1:min(nSamp, size(d,1)), :);

    sig = (d - base) ./ gain;                                 % physical units (mV)
    info = struct('gain', gain, 'baseline', base, 'desc', {desc}, 'nSamp', nSamp);
end
