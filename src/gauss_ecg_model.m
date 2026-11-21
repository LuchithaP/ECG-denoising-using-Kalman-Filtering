function [z, J] = gauss_ecg_model(p, phi)
%GAUSS_ECG_MODEL  Sum-of-Gaussians ECG beat model (McSharry et al., 2003).
%   z(phi) = sum_i alpha_i * exp( -dtheta_i^2 / (2 b_i^2) ),
%   dtheta_i = wrap(phi - theta_i),  i = P, Q, R, S, T.
%   p = [alpha(1:5); b(1:5); theta(1:5)]   (15 x 1)
%   J = dz/dp (numel(phi) x 15) is the analytic Jacobian used by the
%   Levenberg-Marquardt least-squares fit.
    n  = numel(p)/3;
    al = p(1:n)';  b = p(n+1:2*n)';  th = p(2*n+1:3*n)';
    phi = phi(:);
    d  = wrap_pi(phi - th);                 % N x n
    e  = exp(-d.^2 ./ (2*b.^2));
    z  = e * al';
    if nargout > 1
        J = [e, ...                         % d/d alpha
             al .* e .* d.^2 ./ b.^3, ...   % d/d b
             al .* e .* d ./ b.^2];         % d/d theta
    end
end
