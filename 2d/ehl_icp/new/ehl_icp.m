function [R, t, ER, n_iter, Xreg] = ehl_icp(X, Y, opts)
% EHL_ICP: ICP with Hamiltonian rotation update on SO(2)
% and closed-form translation t = q_bar - R x_bar each iteration.
%
% Model used in updates and output:
%   Xreg = R*X + t        (global coordinates)
%
% Inputs:
%   X, Y     : Nx2 source and target
%   opts     : parameters
% Outputs:
%   R, t     : 2x2, 2x1
%   ER       : RMS history (length n_iter+1), ER(1) is initial
%   n_iter   : iterations performed
%   Xreg     : Nx2 registered source


if size(X,2)~=2 || size(Y,2)~=2
    error('X and Y must be Nx2.');
end

% Try original and flipped source; keep best result
candidates = {X};
candidates{end+1} = flipud(X); 

best = []; bestRMS = inf;


for ci = 1:numel(candidates)

    [R1, t1, ER1, n1, Xreg1] = run_once(candidates{ci}, Y, opts);
    
    if ER1(end) < bestRMS
        bestRMS = ER1(end);
        best = struct('R',R1,'t',t1,'ER',ER1,'n',n1,'Xreg',Xreg1);
    end
end

R = best.R; t = best.t; ER = best.ER; n_iter = best.n; Xreg = best.Xreg;

end


% ================== single Hamiltonian run ==================
function [R, t, ER, k, Xreg] = run_once(X, Y, opts)
J = [0 -1; 1  0];

% Initialize global transform
R = eye(2);
t = mean(Y,1)' - R*mean(X,1)';

% Initial RMS on current pose
Xcur = (R*X' + t)';                        % Nx2

idx = knnsearch(Y, Xcur);
d  = Xcur - Y(idx,:);         % Nx2
d2 = sum(d.^2,2);

ER = zeros(opts.max_iter+1,1);
ER(1) = sqrt(mean(d2));

% Momentum state (angular velocity)
v = 0;

for k = 1:opts.max_iter

    % 1) Correspondences on current pose
    Xcur = (R*X' + t)';                    % Nx2
    idx = knnsearch(Y, Xcur);
    Z = Y(idx,:);                          % Nx2

    P = X; 
    Q = Z;

    % 2) Rotation gradient on SO(2) using centered pairs
    x_bar = mean(P,1)';   % 2x1
    q_bar = mean(Q,1)';   % 2x1

    Xc = (P - x_bar')';                       % 2xM
    Qc = (Q - q_bar')';                       % 2xM
    E  = R*Xc - Qc;                           % 2xM residual (centered)
    RJX = R * (J * Xc);                       % 2xM
    theta_grad = (2/size(Xc,2)) * sum(sum(E .* RJX));  % scalar

    % 3) Hamiltonian update (heavy-ball momentum) with clamp
    v = opts.mu * v - opts.eta * theta_grad;          % update "angular momentum"
    step = max(-opts.theta_max, min(opts.theta_max, v));
    R_trial = R * [ cos(step) -sin(step); sin(step) cos(step) ];
    t_trial = q_bar - R_trial * x_bar;                % closed-form t for current matches

    % 4) Optional backtracking on angle step for monotone decrease
    if opts.backtrack
        Xtrial = (R_trial*X' + t_trial)';             % Nx2

        idx = knnsearch(Y, Xtrial);
        d  = Xtrial - Y(idx,:);         % Nx2
        d2trial = sum(d.^2,2);

        f_trial = sqrt(mean(d2trial));
        f_curr  = ER(k);
        bt_cnt = 0;
        while f_trial > f_curr && bt_cnt < opts.bt_max
            step   = step * opts.bt_gamma;
            R_tbt  = R * [ cos(step) -sin(step); sin(step) cos(step) ];
            t_tbt  = q_bar - R_tbt * x_bar;
            Xtrial = (R_tbt*X' + t_tbt)';             % Nx2

            idx = knnsearch(Y, Xtrial);
            d  = Xtrial - Y(idx,:);         % Nx2
            d2trial = sum(d.^2,2);
                        
            f_trial = sqrt(mean(d2trial));
            bt_cnt  = bt_cnt + 1;
            R_trial = R_tbt; t_trial = t_tbt;
        end
    end

    % 5) Commit step
    R = R_trial; t = t_trial;

    % 6) Log RMS on new pose
    Xcur = (R*X' + t)';                            % Nx2

    idx = knnsearch(Y, Xcur);
    d  = Xcur - Y(idx,:);         % Nx2
    d2 = sum(d.^2,2);    
    
    ER(k+1) = sqrt(mean(d2));

    % 7) Stop if converged
%     if abs(ER(k+1)-ER(k)) <= opts.tol * max(1,ER(k))
    if abs((ER(k+1)-ER(k))/ER(k+1)) < opts.tol
        ER = ER(1:k+1);
        break;
    end
end

Xreg = (R*X' + t)';                                % Nx2
end
