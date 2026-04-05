function [result, info] = cg_se3_registration(X1, X2, opts)
% Robust nonmonotone Riemannian CG registration on SE(3) with selectable retraction
% Inputs:
%   X1, X2: Nx3 point clouds (source, target)
%   opts: struct with fields:
%       .max_iter, .xtol, .gtol, .ftol, .alpha_max, .alpha_min, .armijo_eta, .armijo_gamma, .armijo_maxiter
%       .beta_type, .retraction_type ('exp', 'cayley'), .record, etc.
% Outputs:
%   result: struct with .R, .t, .X1_reg, .ER, .k
%   info: struct with .ER, .k
 
N = size(X1,1);
x_c = mean(X1,1)';
y_c = mean(X2,1)';
X1_c = X1 - x_c';   % mean centering
 
% Initialization
R = eye(3);
t = y_c - R*x_c;
ER = zeros(opts.max_iter,1);
k = 0;
 
% Nonmonotone line search params
Q = 1; Cval = Inf; gamma = 0.85;
if isfield(opts, 'gamma'), gamma = opts.gamma; end
 
% ICP initial
X1_transformed = (R * X1_c') + t; X1_transformed = X1_transformed';
% Z = icp_match(X1_transformed, X2);
idx = knnsearch(X2, X1_transformed);    % find closest in target for each point in X1_transformed
Z = X2(idx,:);
 
prev_grad_R = zeros(3,1); prev_grad_t = zeros(3,1);
prev_eta_R = zeros(3,1); prev_eta_t = zeros(3,1);
 
max_rot = 1.0; max_t = 1.0;

while k < opts.max_iter
    k = k+1;
    [err, grad_R, grad_t] = registration_error_se3(R, t, X1_c, Z);
    ER(k) = err;
 
    % --- CG direction update ---
    if k == 1
        eta_R = -grad_R;
        eta_t = -grad_t;
    else
        % Dai-Yuan CG update for SE(3)
        num = grad_R'*grad_R + grad_t'*grad_t;
        den = (grad_R - prev_grad_R)'*prev_eta_R + (grad_t - prev_grad_t)'*prev_eta_t;
        beta = num / max(den, 1e-12);
        if beta < 0 || isnan(beta) || isinf(beta), beta = 0; end
        eta_R = -grad_R + beta * prev_eta_R;
        eta_t = -grad_t + beta * prev_eta_t;
    end
 
    % Clamp directions
    if norm(eta_R) > max_rot, eta_R = eta_R * max_rot / norm(eta_R); end
    if norm(eta_t) > max_t, eta_t = eta_t * max_t / norm(eta_t); end
 
    % --- Nonmonotone Armijo line search ---
    alpha = opts.alpha_max;
    dphi = grad_R'*eta_R + grad_t'*eta_t;
    nls = 1; found = false;
    if k == 1, Cval = err; end
    while nls <= opts.armijo_maxiter
        % Retraction update
        switch lower(opts.retraction_type)
            case 'exp'
                R_new = se3_exponential_retraction(R, eta_R, alpha);
            case 'cayley'
                R_new = se3_cayley_retraction(R, eta_R, alpha);
            otherwise
                error('Unknown retraction_type: %s', opts.retraction_type);
        end
        t_new = t + alpha*eta_t;
 
        X1_transformed_new = (R_new * X1_c') + t_new; X1_transformed_new = X1_transformed_new';
%         Z_new = icp_match(X1_transformed_new, X2);
        idx = knnsearch(X2, X1_transformed_new);    
        Z_new = X2(idx,:);
        err_new = sqrt(mean(sum((X1_transformed_new - Z_new).^2,2)));
 
        % Nonmonotone criterion (GLL style)
        if err_new <= Cval - opts.armijo_eta * alpha * dphi || nls == opts.armijo_maxiter
            found = true; break;
        end
        alpha = alpha * opts.armijo_gamma;
        if alpha < opts.alpha_min, break; end
        nls = nls+1;
    end
 
    % If line search fails, reset to gradient direction
    if ~found || err_new > err
        eta_R = -grad_R; eta_t = -grad_t; alpha = opts.alpha_max * 0.5;
        switch lower(opts.retraction_type)
            case 'exp'
                R_new = se3_exponential_retraction(R, eta_R, alpha);
            case 'cayley'
                R_new = se3_cayley_retraction(R, eta_R, alpha);
        end
        t_new = t + alpha*eta_t;
        X1_transformed_new = (R_new * X1_c') + t_new; X1_transformed_new = X1_transformed_new';
%         Z_new = icp_match(X1_transformed_new, X2);
        idx = knnsearch(X2, X1_transformed_new);    
        Z_new = X2(idx,:);
        err_new = sqrt(mean(sum((X1_transformed_new - Z_new).^2,2)));
    end
 
    % Early termination for NaN/Inf
    if isnan(err_new) || isinf(err_new)
        warning('NaN or Inf encountered, stopping.');
        ER = ER(1:k);
        break;
    end
 
    % --- Stopping criteria ---
    if k > 1
        XDiff = norm(R-R_prev,'fro') + norm(t-t_prev);
        if XDiff < opts.xtol || abs(ER(k)-ER(k-1)) < opts.ftol || err_new < opts.gtol
            ER = ER(1:k);
            break;
        end
    end
 
    % --- Update nonmonotone Q/Cval ---
    Qp = Q; Q = gamma*Qp + 1; Cval = (gamma*Qp*Cval + err_new)/Q;
 
    % --- Store for next iteration ---
    prev_eta_R = eta_R;
    prev_eta_t = eta_t;
    prev_grad_R = grad_R;
    prev_grad_t = grad_t;
    R_prev = R; t_prev = t;
    R = R_new; t = t_new; Z = Z_new;
end
 
result.R = R;
result.t = t;
result.X1_reg = (R * X1' + t)';
result.ER = ER;
result.k = k;
info.ER = ER;
info.k = k;
end



function R_new = se3_cayley_retraction(R, eta_R, alpha)
% Cayley transform retraction for SO(3)

J = [    0     -eta_R(3)  eta_R(2);
      eta_R(3)     0     -eta_R(1);
     -eta_R(2)  eta_R(1)     0   ];

% max_rot = 0.9;
% if norm(eta_R) > max_rot, eta_R = eta_R * max_rot / norm(eta_R); end
epsilon = 1e-8;
Cayley = (eye(3) - 0.5*alpha*J + epsilon*eye(3)) \ (eye(3) + 0.5*alpha*J + epsilon*eye(3));
R_new = R * Cayley;
end


function R_new = se3_exponential_retraction(R, eta_R, alpha)
% Exponential map retraction for SO(3)

J = [    0     -eta_R(3)  eta_R(2);
      eta_R(3)     0     -eta_R(1);
     -eta_R(2)  eta_R(1)     0   ];

% max_rot = 0.2;
% if norm(eta_R) > max_rot, eta_R = eta_R * max_rot / norm(eta_R); end
R_new = R * expm(alpha * J);
end