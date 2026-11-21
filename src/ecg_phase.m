function [phi, omega] = ecg_phase(L, pk, fs)
%ECG_PHASE  Cardiac phase signal from R-peak locations.
%   phi(k) increases linearly from 0 (at an R peak) to 2*pi (at the next R
%   peak) and is wrapped to [-pi, pi).  omega(k) = 2*pi/RR is the local
%   angular frequency in rad/s.  Outside the first/last peaks the nearest
%   RR interval is extrapolated.
    pk = pk(:);
    phi = zeros(L,1); omega = zeros(L,1);
    k = (1:L)';
    % extend peaks at both ends using the neighbouring RR interval
    rr1 = pk(2) - pk(1);  rrN = pk(end) - pk(end-1);
    ext = [flipud((pk(1)-rr1:-rr1:1-rr1)'); pk; (pk(end)+rrN:rrN:L+rrN)'];
    for i = 1:numel(ext)-1
        m = k >= ext(i) & k < ext(i+1);
        rr = ext(i+1) - ext(i);
        phi(m)   = 2*pi*(k(m) - ext(i))/rr;
        omega(m) = 2*pi*fs/rr;
    end
    phi = wrap_pi(phi);
end
