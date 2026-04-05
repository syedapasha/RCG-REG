clear all; close all; clc; 

% shape = "lizard";
% shape = "dog";
% shape = "butterfly";
shape = "bird";


% Load source and target (Nx2)
srcFile = strcat('../../data/',shape,'-s.mat');
tgtFile = strcat('../../data/',shape,'-t.mat');

src = load(srcFile); X1 = src.contour;  % source     
tgt = load(tgtFile); X2 = tgt.contour;  % target     

fprintf('Loaded: source=%d pts, target=%d pts\n', size(X1,1), size(X2,1));

% % Optional but recommended on closed contours: resample to common N
% N = max(size(X1,1), size(X2,1));
% X1 = resample_by_arclength(X1, N, true);
% X2 = resample_by_arclength(X2, N, true);
% fprintf('Resampled to N=%d each\n', N);


% Hamiltonian ICP options
opts.max_iter   = 200;      % max # iterations
opts.tol        = 1e-5;     % stopping threshold

% Hamiltonian (rotation) hyperparameters
opts.eta        = 0.10;     % angle learning rate (start in 0.06–0.12)
opts.mu         = 0.35;     % momentum (0.2–0.5 is typical)
opts.theta_max  = 0.05;     % clamp per iter (rad) ~3°
opts.backtrack  = true;     % Armijo-like backtracking on the angle step
opts.bt_gamma   = 0.7;      % shrink factor for backtracking
opts.bt_max     = 20;       % max backtracking steps

% Try both orientations; keep best
opts.try_flip   = true;


tic
% Register (Hamiltonian-only implementation)
[R, t, ER, k, Xreg] = ehl_icp(X1, X2, opts);
toc


if shape == "lizard"
    lizard_ER_ehl_icp = ER;
    save('lizard_ER_ehl_icp.mat', "lizard_ER_ehl_icp");
elseif shape == "dog"
    dog_ER_ehl_icp = ER;
    save('dog_ER_ehl_icp.mat', "dog_ER_ehl_icp");
elseif shape == "butterfly"
    butterfly_ER_ehl_icp = ER;
    save('butterfly_ER_ehl_icp.mat', "butterfly_ER_ehl_icp");
elseif shape == "bird"
    bird_ER_ehl_icp = ER;
    save('bird_ER_ehl_icp.mat', "bird_ER_ehl_icp");
end


%% --- plot registered shape ---
f = figure('visible', 'off');
plot(X2(:,1), X2(:,2), Xreg(:,1), Xreg(:,2), 'LineStyle','none',  'Marker','.', 'MarkerSize',2);
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
title('EHL-ICP', 'FontSize',15);
legend('target', 'registered', 'FontSize',10, 'Location','southwest');
grid on;
exportgraphics(f, strcat(shape, '_ehl_icp.png'));



fprintf('Iterations: %d\n', k);
fprintf('Final RMS: %.6f\n', ER(end));