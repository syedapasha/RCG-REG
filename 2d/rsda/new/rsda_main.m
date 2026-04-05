clc; clear; close all;

% shape = "lizard";
% shape = "dog";
% shape = "butterfly";
shape = "bird";


% ---------------- Files ----------------
srcFile = strcat('../../data/',shape,'-s.mat');
tgtFile = strcat('../../data/',shape,'-t.mat');

S = load(srcFile); X1 = S.contour;  % source     
T = load(tgtFile); X2 = T.contour;  % target     

fprintf('Loaded: source=%d pts, target=%d pts\n', size(X1,1), size(X2,1));



% ---------------- RSD options ----------------
opts.max_iter       = 300;
% opts.xtol           = 1e-12;
% opts.gtol           = 1e-12;
opts.ftol           = 1e-5; %1e-12;
opts.alpha_max      = 0.3;
opts.alpha_min      = 1e-6;
opts.armijo_eta     = 1e-3;
opts.armijo_gamma   = 0.9;
opts.armijo_maxiter = 20; %12;
opts.retraction_type= 'exp';
opts.theta_max      = 0.03;


% ---------------- Run RSD ----------------
tic
[result, info] = rsd_se2_registration(X1, X2, opts);
toc

Xreg  = result.X1_reg;                         % registered 
theta = atan2(result.R(2,1), result.R(1,1));   % solver rotation (rad)



if shape == "lizard"
    lizard_ER_rsda = result.ER;
    save('lizard_ER_rsda.mat', "lizard_ER_rsda");
elseif shape == "dog"
    dog_ER_rsda = result.ER;
    save('dog_ER_rsda.mat', "dog_ER_rsda");
elseif shape == "butterfly"
    butterfly_ER_rsda = result.ER;
    save('butterfly_ER_rsda.mat', "butterfly_ER_rsda");
elseif shape == "bird"
    bird_ER_rsda = result.ER;
    save('bird_ER_rsda.mat', "bird_ER_rsda");
end


%% --- plot registered shape ---
f = figure('visible', 'off');
plot(X2(:,1), X2(:,2), Xreg(:,1), Xreg(:,2), 'LineStyle','none',  'Marker','.', 'MarkerSize',2);
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
title('RSDA', 'FontSize',15);
legend('target', 'registered', 'FontSize',10, 'Location','southwest');
grid on;
exportgraphics(f, strcat(shape, '_rsda.png'));



% ---------------- Report rotation/translation ----------------
fprintf('\nIterations: %d, Final RMS: %.6f\n', info.k, info.ER(end));
fprintf('Rotation(deg)=%.3f, t=[%.3f %.3f]\n', rad2deg(theta), result.t(1), result.t(2));
