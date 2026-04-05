function [result, info] = rsd_se2_registration(X1, X2, opts)
% Riemannian Steepest Descent (RSD) registration on SE(2) for 2D point sets.
% Matches your RCG solver conventions:
% - Use absolute translation t (not relative to x_c).
% - Transform with Y = R*X1_c + t, X1_c = X1 - x_c.
% - Matching, Armijo backtracking, rotation clamp.
%
% Inputs:
%   X1, X2 : Nx2 point sets (contours), typically closed and resampled to same N.
%
% Outputs:
%   result.R (2x2), result.t (2x1), result.T (3x3), result.X1_reg (Nx2), result.ER (k x 1), result.k
%   info.ER, info.k, info.R, info.t, info.T, info.retraction_type


% ---- checks
assert(size(X1,2)==2 && size(X2,2)==2, 'X1 and X2 must be Nx2.');


% ---- center X1 for rotation; keep X2 absolute
x_c  = mean(X1,1)'; 
X1_c = (X1 - x_c'); 


% ---- initialize (absolute t at target centroid)
R = eye(2);
t = mean(X2,1)';                % absolute translation (identity case stable)

ER = zeros(opts.max_iter,1);
k  = 0;

% ---- initial correspondences (use all points first)
Y = (R * X1_c' + t)';           % absolute coordinates


while k < opts.max_iter
    k = k + 1;

    % Matches 
    idx = knnsearch(X2, Y);
    Z = X2(idx,:);                          

    Xi = X1_c;     % centered source used for gradient
    Zi = Z;                     % matched targets (same rows as Xi)

    % RMS and gradients at (R,t)
    [err, theta_grad, grad_t] = registration_error(R, t, Xi, Zi);
    ER(k) = err;

    % Steepest descent search direction (with simple clamping)
    eta_rot = -theta_grad;
    eta_t   = -grad_t;

    % Normalize directions to avoid huge steps
    if abs(eta_rot) > 1.0, eta_rot = sign(eta_rot); end
    nEt = norm(eta_t); if nEt > 1.0, eta_t = eta_t / nEt; end

    % Descent check: if non-descent, zero out rotation/translation as needed
    dphi = theta_grad*eta_rot + sum(grad_t.*eta_t);
    if dphi >= 0
        % fall back to pure negative gradient (already is), or damp
        eta_rot = -theta_grad; eta_t = -grad_t;
        dphi = theta_grad*eta_rot + sum(grad_t.*eta_t);
        if dphi >= 0
            % If still not descent (degenerate), skip update
            ER = ER(1:k); break;
        end
    end

    % Armijo backtracking with rotation clamp
    alpha = opts.alpha_max;
    if abs(eta_rot) > 0, alpha = min(alpha, opts.theta_max/abs(eta_rot)); end

    accepted = false; err_new = err;
    for ls = 1:opts.armijo_maxiter
        R_new = retract(R, alpha*eta_rot, opts.retraction_type);
        t_new = t + alpha*eta_t;

        Y_trial = (R_new * X1_c' + t_new)';               % absolute coords

        idx = knnsearch(X2, Y_trial);
        Z_trial = X2(idx,:);                          
                
        f1 = sqrt(mean(sum((Y_trial - Z_trial).^2, 2)));

        % Armijo: f(x+a p) <= f(x) + c a dphi (with dphi <= 0)
        if f1 <= err + opts.armijo_eta * alpha * dphi || ls == opts.armijo_maxiter
            accepted = true; err_new = f1;
            break;
        end
        alpha = alpha * opts.armijo_gamma;
        if abs(eta_rot) > 0, alpha = min(alpha, opts.theta_max/abs(eta_rot)); end
        if alpha < opts.alpha_min, break; end
    end

    if ~accepted
        % Very conservative fallback: half-step steepest
        alpha = 0.5 * opts.alpha_max;
        if abs(eta_rot) > 0, alpha = min(alpha, opts.theta_max/abs(eta_rot)); end
        R_new = retract(R, alpha*eta_rot, opts.retraction_type);
        t_new = t + alpha*eta_t;
        Y_trial = (R_new * X1_c' + t_new)'; 

        idx = knnsearch(X2, Y_trial);
        Z_trial = X2(idx,:);                          
        
        err_new = sqrt(mean(sum((Y_trial - Z_trial).^2, 2)));
    end

    % Stopping on small change
    if k > 1
%         XDiff = norm(R-R_new,'fro') + norm(t-t_new);
%         if XDiff < opts.xtol || abs(ER(k)-ER(k-1)) < opts.ftol || err_new < opts.gtol
        if abs((ER(k)-ER(k-1))/ER(k)) < opts.ftol 
            R = R_new; t = t_new; ER = ER(1:k); break;
        end
    end

    % Commit and refresh
    R = R_new; t = t_new;
    Y = (R * X1_c' + t)';       % refresh transformed points
end

% ---- outputs
X1_reg = (R * X1_c' + t)';      % absolute coords
result.R = R; result.t = t; result.T = [R, t; 0 0 1];
result.X1_reg = X1_reg; result.X_reg = X1_reg;
result.ER = ER(1:k); result.k = k;

info.ER = result.ER; info.k = k; info.R = R; info.t = t; info.T = result.T;
info.retraction_type = opts.retraction_type;



% ====================== helpers ======================


% ***********************************************
% retraction on SE(2).
%
% Inputs:
%   R      : rotation matrix 
%   dtheta : small angle rotation
%   how    : retraction mappng (e.g. exp, Cayley).
%
% Outputs:
%   Rn     : retraction 
% ***********************************************
    
function Rn = retract(R, dtheta, how)
switch lower(how)
    case 'exp'
        Rn = R * [ cos(dtheta) -sin(dtheta); sin(dtheta) cos(dtheta) ];
    case 'cayley'
        ang = 2*atan(0.5*dtheta);
        Rn = R * [ cos(ang) -sin(ang); sin(ang) cos(ang) ];
    otherwise
        Rn = R * [ cos(dtheta) -sin(dtheta); sin(dtheta) cos(dtheta) ];
end
end


% ***************************************************************
% registration error for 2D point sets.
%
% Inputs:
%   R, t : rotation matrix and translation vector
%   Xi, Zi : Nx2 point sets (contours).
%
% Outputs:
%   f       : rms error 
%   gth, gt : gradients wrt rotation and translation respectively
% ***************************************************************

function [f, gth, gt] = registration_error(R, t, Xi, Zi)

J = [0 -1; 1 0];
gt = [0; 0]; gth = 0; acc = 0;

Nm = size(Xi,1);
for i = 1:Nm
    xi = Xi(i,:)'; zi = Zi(i,:)';
    ei = R*xi + t - zi;
    acc = acc + ei.'*ei;
    gt  = gt  + ei;
    gth = gth + ei.'*(R*(J*xi));
end

f   = sqrt(acc / Nm);
gt  = (2/Nm) * gt;
gth = (2/Nm) * gth;
end


end