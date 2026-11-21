function [p, hist] = fit_gaussians_lm(phc, mu, maxIter)
%FIT_GAUSSIANS_LM  Nonlinear least-squares fit of the 5-Gaussian ECG model.
%   Minimises  J(p) = sum_j ( mu(phi_j) - z(phi_j; p) )^2  with the
%   Levenberg-Marquardt algorithm (implemented here, no toolbox needed):
%       (J'J + lambda*diag(J'J)) dp = J' r
%   lambda is decreased after a successful step (Gauss-Newton like) and
%   increased after a failed one (gradient-descent like).
%
%   phc, mu : phase bins and mean beat (from phase_mean_beat)
%   p       : [alpha; b; theta] for waves P, Q, R, S, T
%   hist    : struct with cost and parameter history (for plots)

    if nargin < 3, maxIter = 100; end
    phc = phc(:); mu = mu(:);

    % ---------- data-driven initial guess ----------
    pick = @(lo, hi, f) local_ext(phc, mu, lo, hi, f);
    [thR, aR] = pick(-0.15, 0.15, @(v) max(abs(v)));   % R wave (either polarity)
    pol = sign(aR);  if pol == 0, pol = 1; end
    neg = @(v) max(-pol*v);                             % deflection opposite to R
    [thQ, aQ] = pick(thR-0.35, thR-0.04, neg);  aQ = -pol*aQ;
    [thS, aS] = pick(thR+0.04, thR+0.35, neg);  aS = -pol*aS;
    [thP, aP] = pick(-2.6, thR-0.5, @max);
    [thT, aT] = pick(thR+0.6, 2.9, @max);
    al0 = [aP; aQ; aR; aS; aT];
    b0  = [0.20; 0.05; 0.05; 0.05; 0.35];
    th0 = [thP; thQ; thR; thS; thT];

    % Physiological box constraints (projected Levenberg-Marquardt).
    % Without them the fit can be "mathematically good but physiologically
    % wrong" (e.g. huge R and S waves cancelling each other).
    %   order   : P < Q < R < S < T around the R peak
    %   widths  : narrow Q, R, S (QRS complex), broader P and T
    %   signs   : Q and S deflect opposite to R
    bLo  = [0.08; 0.02; 0.02; 0.02; 0.15];
    bHi  = [0.50; 0.12; 0.12; 0.12; 0.90];
    thLo = [-2.9; thR-0.35; thR-0.05; thR+0.04; thR+0.5];
    thHi = [thR-0.4; thR-0.04; thR+0.05; thR+0.35; 3.0];
    sgn  = @(a) [a(1); -pol*max(-pol*a(2),0); pol*max(pol*a(3),0); -pol*max(-pol*a(4),0); a(5)];
    project = @(q) [sgn(q(1:5)); min(max(q(6:10), bLo), bHi); min(max(q(11:15), thLo), thHi)];
    p = project([al0; b0; th0]);

    lambda = 1e-2;
    r = mu - gauss_ecg_model(p, phc);
    cost = r'*r;
    hist.cost = cost; hist.p = p;
    for it = 1:maxIter
        [z, J] = gauss_ecg_model(p, phc);
        r  = mu - z;
        H  = J'*J;  g = J'*r;
        improved = false;
        while ~improved && lambda < 1e10
            % small ridge term keeps the system solvable when a wave amplitude
            % sits at 0 (its width/centre then have no influence on the cost)
            dp   = (H + lambda*diag(diag(H) + 1e-6*trace(H)/15)) \ g;
            pNew = project(p + dp);                        % keep inside the box
            rNew = mu - gauss_ecg_model(pNew, phc);
            cNew = rNew'*rNew;
            if all(isfinite(dp)) && cNew < cost
                p = pNew; lambda = max(lambda/10, 1e-7); improved = true;
                relChange = (cost - cNew)/cost;  cost = cNew;
            else
                lambda = lambda*10;
            end
        end
        hist.cost(end+1) = cost; hist.p(:,end+1) = p;
        if ~improved || relChange < 1e-9, break; end
    end
    hist.p0 = project([al0; b0; th0]);
end

function [th, a] = local_ext(phc, mu, lo, hi, f)
    m = find(phc >= lo & phc <= hi);
    [a, j] = f(mu(m));
    th = phc(m(j));
end
