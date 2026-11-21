function r = ekf_pipeline(xRaw, fs, opt)
%EKF_PIPELINE  Complete model-based denoising chain for one noisy ECG.
%   1) baseline removal (median filters)
%   2) R-peak detection on the NOISY signal -> cardiac phase phi(t)
%   3) LS:  phase-domain mean beat, Levenberg-Marquardt fit of 5 Gaussians
%   4) MLE: measurement noise variance  sigma^2 = mean( (s - z(phi))^2 )
%           (the ML estimate for i.i.d. Gaussian noise)
%   5) EKF forward pass + RTS smoother (EKS)
    if nargin < 3, opt = struct; end
    r.xb  = remove_baseline(xRaw, fs);
    r.pk  = detect_rpeaks(r.xb, fs);
    [r.phi, r.omega] = ecg_phase(numel(r.xb), r.pk, fs);
    [r.phc, r.mu, r.sd] = phase_mean_beat(r.xb, r.phi, 250);
    [r.p, r.lsHist] = fit_gaussians_lm(r.phc, r.mu);

    % ML estimate of the ECG measurement noise variance
    resid    = r.xb - gauss_ecg_model(r.p, r.phi);
    r.sigma2 = mean(resid.^2);
    sigPhi2  = (2*pi*0.02)^2;                 % phase observation noise
    R = diag([sigPhi2, r.sigma2]);

    r.ekf = ekf_ecg(r.xb, r.phi, r.omega, r.p, fs, R, opt);
    r.xEKF = r.ekf.zf;
    r.xEKS = r.ekf.zs;
end
