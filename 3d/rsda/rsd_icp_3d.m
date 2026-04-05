function [R, t, errors, n_iter, Xreg] = rsd_icp_3d(X, Y, varargin)
opts = inputParser;
addParameter(opts, 'max_iter', 30);
addParameter(opts, 'min_iter', 3);
addParameter(opts, 'eta', 01);
addParameter(opts, 'tol', 1e-5);
parse(opts, varargin{:});
max_iter = opts.Results.max_iter;
min_iter = opts.Results.min_iter;
eta = opts.Results.eta;
tol = opts.Results.tol;

d = size(X,1); N = size(X,2);
if d ~= 3, error('X and Y must be 3xN.'); end
if size(Y,1) ~= d, error('X and Y must have the same dimensionality.'); end

% Initial rotation and translation
R = eye(3);
t = mean(Y,2) - mean(X,2);

% Centralize data
x_c = mean(X,2);
y_c = mean(Y,2);

errors = zeros(max_iter+1,1);

% Initial registration and error
Xreg = R*X + t;
[~, corrY] = min(pdist2(Xreg', Y'), [], 2);
err = sqrt(mean(sum((Xreg - Y(:,corrY)).^2, 1)));
errors(1) = err;

for k = 1:max_iter
    k
    % Step 1: Find closest points (correspondence)
    Xreg = R*X + t;
    [~, corrY] = min(pdist2(Xreg', Y'), [], 2);
    Z = Y(:, corrY);

    % Step 2: Compute RMS error
    err = sqrt(mean(sum((Xreg - Z).^2, 1)));
    errors(k+1) = err;

    % Step 3: Centralize
    x_c = mean(X,2);
    z_c = mean(Z,2);
    Xc = X - x_c;
    Zc = Z - z_c;

    % Step 4: Compute gradient wrt rotation (Eq. 26)
    G = zeros(3);
    for i = 1:N
        G = G + Xc(:,i) * (R*Xc(:,i) + t - Zc(:,i))';
    end
    G = (2/N) * G;

    % Step 5: Riemannian gradient step (NO momentum)
    J = R' * G;
    J = 0.5*(J - J'); % Project to skew-symmetric
    R = R * expm(eta * J);

    % Step 6: Update translation
    t = y_c - R * x_c;

    % Step 7: Check convergence (avoid early exit)
    if k >= min_iter
        rel_change = abs(errors(k+1)-errors(k))/max(errors(k),1e-12);
        if rel_change < tol
            break;
        end
    end
end

n_iter = k;
errors = errors(1:n_iter+1);
Xreg = R*X + t;
end