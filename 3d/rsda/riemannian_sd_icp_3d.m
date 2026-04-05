function [R, t, errors, n_iter, Xreg] = riemannian_sd_icp_3d(X1, X2, varargin)
% Riemannian Steepest Descent ICP for 3D rigid registration
%   Inputs:
%     X1, X2: 3 x N source and target point clouds
%     varargin: 'max_iter', 'min_iter', 'eta', 'tol'
%   Outputs:
%     R: 3x3 rotation
%     t: 3x1 translation
%     errors: RMS error history
%     n_iter: number of iterations
%     Xreg: registered source points

p = inputParser;
addParameter(p, 'max_iter', 1000);
addParameter(p, 'min_iter', 3);
addParameter(p, 'eta', 0.2);
addParameter(p, 'tol', 1e-12);
parse(p, varargin{:});

max_iter = p.Results.max_iter;
min_iter = p.Results.min_iter;
eta = p.Results.eta;
tol = p.Results.tol;

N = size(X1,2);

% Initialization
R = eye(3);
t = zeros(3,1);

errors = zeros(1, max_iter+1);

for k = 1:max_iter
    % Transform source
    Xreg = R*X1 + t;

    % ICP: Find closest points in target (nearest neighbor)
    idx = knnsearch(X2', Xreg');
    X2_corr = X2(:,idx);

    % Compute RMS error
    err = sqrt(mean(sum((X2_corr - Xreg).^2,1)));
    errors(k+1) = err;

    % Compute gradients
    dX = Xreg - X2_corr; % 3xN

    % Translation gradient: mean error direction
    grad_t = mean(dX,2);

    % Rotation gradient: sum of cross products
    grad_R = zeros(3,3);
    for i = 1:N
        grad_R = grad_R + dX(:,i) * X1(:,i)';
    end
    grad_R = grad_R / N;

    % Map to so(3) (skew-symmetric part)
    Omega = grad_R - grad_R';
    omega = [Omega(3,2); Omega(1,3); Omega(2,1)]; % Axis-angle

    % Update translation (Euclidean step)
    t = t - eta * grad_t;

    % Update rotation (Exponential map on SO(3))
    theta = norm(omega);
    if theta > 1e-12
        Omega_hat = skew(omega/theta);
        dR = expm(-eta * theta * Omega_hat);
        R = dR * R; % Left-invariant update
    end

    % Convergence check
    if k > min_iter && abs(errors(k+1)-errors(k)) < tol
        break;
    end
end

n_iter = k;
errors = errors(1:n_iter+1);
Xreg = R*X1 + t;

end

function S = skew(w)
% Convert 3x1 vector to 3x3 skew-symmetric matrix
S = [  0   -w(3)  w(2);
      w(3)   0   -w(1);
     -w(2) w(1)    0 ];
end