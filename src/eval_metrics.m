function m = eval_metrics(xClean, xNoisy, xHat, fs)
%EVAL_METRICS  Quantitative comparison of an estimate with the clean ECG.
%   The first and last second are excluded (filter start-up transients).
%   snrIn   : 10 log10( sum x^2 / sum (y - x)^2 )      
%   snrOut  : 10 log10( sum x^2 / sum (xhat - x)^2 )    
%   snrImp  : snrOut - snrIn                            
%   rmse    : root-mean-square error (mV)
%   prd     : percentage root-mean-square difference (%)
%   corr    : Pearson correlation between estimate and clean ECG
    k = (fs+1) : (numel(xClean)-fs);
    x = xClean(k); y = xNoisy(k); e = xHat(k) - x;
    m.snrIn  = 10*log10(sum(x.^2) / sum((y - x).^2));
    m.snrOut = 10*log10(sum(x.^2) / sum(e.^2));
    m.snrImp = m.snrOut - m.snrIn;
    m.rmse   = sqrt(mean(e.^2));
    m.prd    = 100*sqrt(sum(e.^2)/sum(x.^2));
    c = corrcoef(x, xHat(k));  m.corr = c(1,2);
end
