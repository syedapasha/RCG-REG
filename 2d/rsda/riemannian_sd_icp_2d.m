function [R, t, errors, n_iter, Xreg] = riemannian_sd_icp_2d(X1, X2, varargin)
% Riemannian steepest descent ICP for 2D rigid registration (SE(2))

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

R = eye(2);
t = zeros(2,1);

errors = zeros(1, max_iter+1);

% --- compute initial RMS error ---
Xreg = R*X1 + t;                % initial registration
idx = knnsearch(X2', Xreg');    % find closest in target for each point in Xreg
X2_corr = X2(:,idx);
errors(1) = sqrt(mean(sum((X2_corr - Xreg).^2,1)));     % initial RMS error

for k = 1:max_iter
    % ICP iteration
    Xreg = R*X1 + t;
    idx = knnsearch(X2', Xreg');
    X2_corr = X2(:,idx);
    err = sqrt(mean(sum((X2_corr - Xreg).^2,1)));
    errors(k+1) = err;

    % --- gradient descent step ---
    dX = Xreg - X2_corr;    % 2xN
    grad_t = mean(dX,2);

    grad_R = zeros(2,2);
    for i = 1:N
        grad_R = grad_R + dX(:,i) * X1(:,i)';
    end
    grad_R = grad_R / N;

    Omega = grad_R - grad_R';
    omega = Omega(2,1);     % scalar in 2D

    t = t - eta * grad_t;

    theta = abs(omega);
    if theta > 1e-12
        Omega_hat = [0 -1; 1 0];
        dR = expm(-eta * theta * Omega_hat * sign(omega));
        R = dR * R;
    end

    if k > min_iter && abs(errors(k+1)-errors(k)) < tol
        break;
    end
end

n_iter = k;
errors = errors(1:n_iter+1);
Xreg = R*X1 + t;

end