function [samp, type] = read_wfdb_ann(annFile)
%READ_WFDB_ANN  Read a PhysioNet MIT-format annotation file (e.g. 100.atr).
%   [samp, type] = read_wfdb_ann('data/raw/100.atr')
%   samp : annotation sample indices (1-based)
%   type : annotation codes (1 = N normal beat, 5 = V PVC, ... see WFDB docs)

    fid = fopen(annFile, 'r');
    w = fread(fid, inf, 'uint16', 0, 'l');
    fclose(fid);

    samp = zeros(numel(w),1); type = zeros(numel(w),1);
    n = 0; t = 0; i = 1;
    while i <= numel(w)
        A = bitshift(w(i), -10);          % 6-bit annotation code
        I = bitand(w(i), 1023);           % 10-bit interval / length field
        if A == 0 && I == 0, break; end   % end of file
        switch A
            case 59                       % SKIP: 32-bit interval follows
                t = t + bitshift(w(i+1),16) + w(i+2);
                i = i + 3; continue;
            case 63                       % AUX string: skip its bytes
                i = i + 1 + ceil(I/2); continue;
            case {60, 61, 62}             % NUM / SUB / CHAN modifiers
                i = i + 1; continue;
            otherwise                     % ordinary annotation
                t = t + I;
                n = n + 1; samp(n) = t + 1; type(n) = A;
                i = i + 1;
        end
    end
    samp = samp(1:n); type = type(1:n);
end
