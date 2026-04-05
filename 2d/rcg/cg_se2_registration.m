function [result, info] = cg_se2_registration(X1, X2, opts)
% Riemannian Conjugate Gradient (RCG) registration on SE(2) with selectable retraction

% center X1 for rotation; keep X2 in absolute coordinates
x_c = mean(X1,1)';
y_c = mean(X2,1)';
X1_c = X1 - x_c';   % mean centering


%% --- Initialization ---

% initial rotation and translation
R = eye(2);
t = y_c;

% --- initial correspondences ---
X1_transformed = (R * X1_c') + t; X1_transformed = X1_transformed';
% find closest in X_2 for each point in X1_transformed
idx = knnsearch(X2, X1_transformed);    
Z = X2(idx,:);

ER = zeros(opts.max_iter,1);
k = 0;

% For CG
prev_theta_grad = 0; prev_grad_t = [0;0];
prev_eta_rot = 0; prev_eta_t = [0;0];

% max_rot = 0.0087; max_t = 1.0;
max_rot = 0.03; max_t = 1.0;


%% --- Riemannian CG iterations ---

while k < opts.max_iter

    k = k+1;
    
    % RMS + gradients at (R,t): ei = R*xi + t - zi
    [err, theta_grad, grad_t] = registration_error(R, t, X1_c, Z);
    ER(k) = err;

    % --- CG direction update ---
    if k == 1       % gradient descent direction
        eta_rot = -theta_grad;
        eta_t   = -grad_t;
    else            % Dai–Yuan CG direction 
        num = theta_grad^2 + sum(grad_t.^2);
        den = (theta_grad - prev_theta_grad)*prev_eta_rot ...
            + sum((grad_t-prev_grad_t).*prev_eta_t);
        
        if abs(den) < 1e-12 || isnan(den) || isinf(den)
            beta = 0;
        else
            beta = num / den;
            if beta < 0 || isnan(beta) || isinf(beta), beta = 0; end
        end

        eta_rot = -theta_grad + beta * prev_eta_rot;
        eta_t   = -grad_t     + beta * prev_eta_t;
    end

    % Clamp direction to avoid wild updates
    if abs(eta_rot) > max_rot, eta_rot = sign(eta_rot) * max_rot; end
    if norm(eta_t) > max_t, eta_t = eta_t * max_t / norm(eta_t); end


    % --- Armijo backtracking line search ---
    alpha = opts.alpha_max;

    for ls = 1:opts.armijo_maxiter

        % select retraction type for rotation
        switch lower(opts.retraction_type)
            case 'exp'
                R_new = se2_exponential_retraction(R, eta_rot, alpha);
            case 'cayley'
                R_new = se2_cayley_retraction(R, eta_rot, alpha);
            case 'identity'
                R_new = se2_identity_retraction(R, eta_rot, alpha);
            otherwise
                error('Unknown retraction_type: %s', opts.retraction_type);
        end
        t_new = t + alpha * eta_t;

        X1_transformed_new = (R_new * X1_c') + t_new; X1_transformed_new = X1_transformed_new';
        idx = knnsearch(X2, X1_transformed_new);    
        Z_new = X2(idx,:);
        err_new = sqrt(mean(sum((X1_transformed_new - Z_new).^2,2)));
        
        % Armijo condition
        dphi = theta_grad*eta_rot + sum(grad_t.*eta_t);
        if err_new <= err - opts.armijo_eta * alpha * dphi || ls == opts.armijo_maxiter
            break;
        end
        alpha = alpha * opts.armijo_gamma;
        if alpha < opts.alpha_min, break; end
    end


    % If line search fails (error increases), reset direction to gradient
    if err_new > err
        eta_rot = -theta_grad;
        eta_t   = -grad_t;

        % Clamp direction to avoid wild updates
        if abs(eta_rot) > max_rot, eta_rot = sign(eta_rot) * max_rot; end
%         if norm(eta_t) > max_t, eta_t = eta_t * max_t / norm(eta_t); end
        
        alpha = opts.alpha_max * 0.5;
        
        switch lower(opts.retraction_type)
            case 'exp'
                R_new = se2_exponential_retraction(R, eta_rot, alpha);
            case 'cayley'
                R_new = se2_cayley_retraction(R, eta_rot, alpha);
            case 'identity'
                R_new = se2_identity_retraction(R, eta_rot, alpha);
        end

        t_new = t + alpha * eta_t;

        X1_transformed_new = (R_new * X1_c') + t_new; X1_transformed_new = X1_transformed_new';
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

    % --- Store for next iteration ---
    prev_eta_rot = eta_rot;
    prev_eta_t = eta_t;
    prev_theta_grad = theta_grad;
    prev_grad_t = grad_t;
    R_prev = R; t_prev = t;
    R = R_new; t = t_new; Z = Z_new;
end

result.R = R;
result.t = t;
result.Xreg = (R * X1_c' + t)';
result.ER = ER;
result.k = k;
info.ER = ER;
info.k = k;
end




%% --- Exponential Retraction ---
function R_new = se2_exponential_retraction(R, eta_rot, alpha)
% Exponential map retraction for SO(2) with clamped step size

cos_th = cos(alpha*eta_rot);
sin_th = sin(alpha*eta_rot);

R_new = R * [ cos_th, -sin_th; 
              sin_th,  cos_th ];
end


%% --- Cayley Retraction ---
function R_new = se2_cayley_retraction(R, eta_rot, alpha)
% Cayley transform retraction for SO(2)
J = [0 -1; 1 0];

epsilon = 1e-6; % Regularization for stability
Cayley = (eye(2) - 0.5*alpha*eta_rot*J + epsilon*eye(2)) \ (eye(2) + 0.5*alpha*eta_rot*J + epsilon*eye(2));
R_new = R * Cayley;
end


%% --- Identity Retraction ---
function R_new = se2_identity_retraction(R, eta_rot, alpha)
% Identity retraction for SO(2) with clamped step size and projection to SO(2)
J = [0 -1; 1 0];

R_new = R * (eye(2) + alpha*eta_rot*J);
[U,~,V] = svd(R_new);
R_new = U*V'; % Project back to SO(2)
end