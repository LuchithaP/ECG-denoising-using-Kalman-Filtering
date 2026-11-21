function [phc, mu, sd] = phase_mean_beat(x, phi, nBins)
%PHASE_MEAN_BEAT  Ensemble-average beat in the phase domain.
%   Every sample is assigned to one of nBins phase bins in [-pi, pi);
%   mu and sd are the mean and standard deviation of the ECG in each bin.
%   Averaging over many beats suppresses noise by ~sqrt(#beats).
    if nargin < 3, nBins = 250; end
    edges = linspace(-pi, pi, nBins+1);
    phc   = (edges(1:end-1) + edges(2:end))'/2;
    bin   = discretize(phi, edges);
    mu = accumarray(bin, x, [nBins 1], @mean, NaN);
    sd = accumarray(bin, x, [nBins 1], @std,  NaN);
    ok = ~isnan(mu);                                % fill any empty bins
    mu = interp1(phc(ok), mu(ok), phc, 'linear', 'extrap');
    sd = interp1(phc(ok), sd(ok), phc, 'linear', 'extrap');
end
