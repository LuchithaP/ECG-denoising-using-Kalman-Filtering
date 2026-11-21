function v = make_noise(ds, type, snrDb, L, offset, P)
%MAKE_NOISE  Real NSTDB noise scaled to a target nst-SNR.
%   v = make_noise(ds, type, snrDb, L, offset, P)
%   type   : 'ma' (muscle artifact), 'bw' (baseline wander) or 'ma+bw'
%   snrDb  : target SNR in dB (PhysioNet nst definition)
%   L      : number of samples;  offset : start sample in the noise record
%   P      : struct with fields S, N_ma, N_bw (from nst_power)
%
%   Gain:  SNR = 10 log10( S / (N a^2) )  ->  a = sqrt( S / (N 10^(SNR/10)) )
%   For 'ma+bw' the total noise power budget S/10^(SNR/10) is split equally
%   between the two sources (power-based allocation), i.e. each gets half.

    seg = @(x) x(offset:offset+L-1) - mean(x(offset:offset+L-1));
    Ntarget = P.S / 10^(snrDb/10);                  % allowed noise power
    switch lower(type)
        case 'ma'
            v = sqrt(Ntarget / P.N_ma) * seg(ds.noise_ma);
        case 'bw'
            v = sqrt(Ntarget / P.N_bw) * seg(ds.noise_bw);
        case 'ma+bw'
            v = sqrt(0.5*Ntarget / P.N_ma) * seg(ds.noise_ma) + ...
                sqrt(0.5*Ntarget / P.N_bw) * seg(ds.noise_bw);
        otherwise
            error('Unknown noise type %s', type);
    end
end
