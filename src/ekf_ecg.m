function out = ekf_ecg(s, phi, omega, p, fs, R, opt)
%EKF_ECG  Extended Kalman Filter + RTS smoother for ECG denoising.
%   Implements the nonlinear Bayesian ECG filter of Sameni et al. (2007),
%   built on the McSharry dynamical ECG model.
%
%   STATE      x_k = [theta_k ; z_k]   (cardiac phase, clean ECG amplitude)
%   PROCESS    theta_{k+1} = wrap( theta_k + omega*dt )
%              z_{k+1}     = z_k - sum_i dt*alpha_i*omega/b_i^2 * dth_i
%                                 * exp(-dth_i^2/(2 b_i^2)) + eta
%              dth_i = wrap(theta_k - theta_i)
%              Process noise w = [alpha(5) b(5) theta(5) omega eta] (17 x 1):
%              the wave parameters are random around their LS estimates,
%              which lets the model follow beat-to-beat changes.
%   MEASURE    y_k = [phi_k ; s_k] = x_k + v_k ,  v_k ~ N(0, R)
%              phi_k from R peaks, s_k = noisy ECG sample.
%
%   s, phi, omega : noisy ECG, observed phase, angular frequency (rad/s)
%   p             : LS-fitted Gaussian parameters [alpha; b; theta]
%   R             : 2x2 measurement noise covariance
%   opt           : struct with process-noise scales (see defaults below)
%   out.zf / out.zs : EKF (filtered) and EKS (smoothed) ECG estimates

    if nargin < 7, opt = struct; end
    % Process-noise standard deviations (relative to the LS parameters).
    % Tuned by a small grid search on a noise segment NOT used for evaluation.
    % (These defaults equal the result of tune_parameters.m.)
    d = struct('qa', 0.05, 'qb', 0.05, 'qth', 0.02, 'qw', 0.05, 'qeta', 0.001);
    fn = fieldnames(d);
    for i = 1:numel(fn), if ~isfield(opt, fn{i}), opt.(fn{i}) = d.(fn{i}); end, end

    N  = numel(s);  dt = 1/fs;  n = numel(p)/3;
    al = p(1:n);  b = p(n+1:2*n);  th = p(2*n+1:3*n);
    w0 = median(omega);
    Q  = diag([(opt.qa*abs(al)).^2; (opt.qb*b).^2; opt.qth^2*ones(n,1); ...
               (opt.qw*w0)^2; (opt.qeta*max(abs(al)))^2]);

    % storage for the smoother
    xp = zeros(2,N); Pp = zeros(2,2,N);     % predicted
    xf = zeros(2,N); Pf = zeros(2,2,N);     % filtered
    A21 = zeros(1,N);                       % Jacobian dz/dtheta used for step k-1 -> k
    Kz  = zeros(1,N);                       % Kalman gain on the ECG sample

    x = [phi(1); s(1)];  P = diag([0.1, 1]);
    for k = 1:N
        % ---------------- PREDICT ----------------
        if k > 1
            om  = omega(k);
            dth = wrap_pi(x(1) - th);
            e   = exp(-dth.^2 ./ (2*b.^2));
            % state propagation with noise at its mean
            xpk = [wrap_pi(x(1) + om*dt); ...
                   x(2) - dt*sum(al.*om./b.^2 .* dth .* e)];
            % Jacobian wrt state
            a21 = -dt*sum(al.*om./b.^2 .* (1 - dth.^2./b.^2) .* e);
            A   = [1 0; a21 1];
            % Jacobian wrt process noise w
            F = zeros(2, 3*n+2);
            F(1, 3*n+1) = dt;
            F(2, 1:n)       = -dt*om*dth.*e ./ b.^2;
            F(2, n+1:2*n)   =  2*dt*al.*om.*dth.*e ./ b.^3 .* (1 - dth.^2./(2*b.^2));
            F(2, 2*n+1:3*n) =  dt*al.*om ./ b.^2 .* (1 - dth.^2./b.^2) .* e;
            F(2, 3*n+1)     = -dt*sum(al.*dth.*e ./ b.^2);
            F(2, 3*n+2)     =  1;
            Ppk = A*P*A' + F*Q*F';
            A21(k) = a21;
        else
            xpk = x;  Ppk = P;
        end
        xp(:,k) = xpk;  Pp(:,:,k) = Ppk;

        % ---------------- UPDATE -----------------
        nu = [wrap_pi(phi(k) - xpk(1)); s(k) - xpk(2)];    % innovation
        Sk = Ppk + R;
        K  = Ppk / Sk;                                     % Kalman gain
        x  = xpk + K*nu;   x(1) = wrap_pi(x(1));
        P  = (eye(2) - K)*Ppk;
        xf(:,k) = x;  Pf(:,:,k) = P;  Kz(k) = K(2,2);
    end

    % ------------- RTS SMOOTHER (backward pass) -------------
    xs = xf;  Ps = Pf;
    for k = N-1:-1:1
        A  = [1 0; A21(k+1) 1];
        C  = Pf(:,:,k) * A' / Pp(:,:,k+1);
        dx = xs(:,k+1) - xp(:,k+1);  dx(1) = wrap_pi(dx(1));
        xs(:,k) = xf(:,k) + C*dx;    xs(1,k) = wrap_pi(xs(1,k));
        Ps(:,:,k) = Pf(:,:,k) + C*(Ps(:,:,k+1) - Pp(:,:,k+1))*C';
    end

    out.zf = xf(2,:)';  out.zs = xs(2,:)';
    out.thf = xf(1,:)'; out.ths = xs(1,:)';
    out.Pzf = squeeze(Pf(2,2,:));  out.Pzs = squeeze(Ps(2,2,:));
    out.Kz = Kz(:);  out.Q = Q;  out.R = R;
end
