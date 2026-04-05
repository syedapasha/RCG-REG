% Runner for RCG rigid registration on 2D contours (Nx2).
% Uses all points; resamples both to a common N via arclength for stability.

clc; clear; close all;

srcFile = 'lizzard-11.mat';    % source
tgtFile = 'lizzard-13.mat';   % target

S = load(srcFile); X1 = S.contour;
T = load(tgtFile); X2 = T.contour;

fprintf('Loaded: source=%d pts, target=%d pts\n', size(X1,1), size(X2,1));

% Make both have identical N (recommended for closed contours)
N = max(size(X1,1), size(X2,1));
X1 = resample_by_arclength(X1, N, true);
X2 = resample_by_arclength(X2, N, true);
fprintf('Resampled to common N=%d\n', N);

% RCG options (safe for 1500+ points; smoother RMS)
opts.max_iter       = 200;
opts.xtol           = 1e-12;
opts.gtol           = 1e-12;
opts.ftol           = 1e-12;
opts.alpha_max      = 0.3;     % more conservative step
opts.alpha_min      = 1e-6;
opts.beta_type      = 'DY';
opts.armijo_eta     = 1e-3;
opts.armijo_gamma   = 0.9;
opts.armijo_maxiter = 12;
opts.retraction_type= 'exp';
opts.theta_max      = 0.03;    % ~1.7 deg per iter; helps smooth increments

% IMPORTANT for smooth ER curve: avoid changing pair sets during solve
opts.trim_keep      = 1.0;     % no trimming -> fixed-size sum
opts.mutual         = false;   % no mutual-NN gating for smoother logging

[result, info] = cg_se2_registration(X1, X2, opts);

% Plot overlays and RMS curve (requires plot_registration.m and local_nn.m)
plot_registration(X1, X2, result, info, ...
    'ShowNNRMS', true, ...     % prints final NN RMS with fixed mapping
    'SmoothER',  true, ...
    'SmoothAlpha', 0.15, ...
    'ShowCorr',  false);

% Extra: print a fixed-mapping NN RMS (no rematching across iterations)
Xa = result.X1_reg;
try
    [idx_eval, ~] = local_nn(Xa, X2);
    Z_eval = X2(idx_eval, :);
    r_fixed = sqrt(mean(sum((Xa - Z_eval).^2, 2)));
    fprintf('Final fixed-mapping NN RMS (no trim/mutual): %.6f\n', r_fixed);
catch
    warning('local_nn.m not found; fixed-mapping RMS not printed.');
end

fprintf('\nIterations: %d, Line-search RMS (last stored): %.6f\n', numel(info.ER), info.ER(end));