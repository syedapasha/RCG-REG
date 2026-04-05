clear all; close all;

shape = "lizard";
% shape = "dog";
% shape = "butterfly";
% shape = "bird";

%% Load source and target shapes
%*******************************
srcFile = strcat('../data/',shape,'-s.mat');
tgtFile = strcat('../data/',shape,'-t.mat');

S = load(srcFile); X1 = S.contour;  % source
T = load(tgtFile); X2 = T.contour;  % target



opts.max_iter = 40;
opts.retraction_type = 'identity';
opts.xtol = 1e-5;
opts.gtol = 1e-5;
opts.ftol = 1e-5;
opts.alpha_max = 0.9;   
opts.alpha_min = 1e-6;
opts.beta_type = 'DY';   
opts.armijo_eta = 0.5;  
opts.armijo_gamma = 0.9; 
opts.armijo_maxiter = 12;

tic
[result, info] = cg_se2_registration(X1, X2, opts);
toc

Xreg = result.Xreg;
ER = result.ER;
k = result.k;

% lizard_ER_id = ER;
% save('lizard_er_id.mat', "lizard_ER_id");
% return

if shape == "lizard"
    lizard_ER_rcg = ER;
    save('lizard_er_rcg.mat', "lizard_ER_rcg");
elseif shape == "dog"
    dog_ER_rcg = ER;
    save('dog_er_rcg.mat', "dog_ER_rcg");
elseif shape == "butterfly"
    butterfly_ER_rcg = ER;
    save('butterfly_er_rcg.mat', "butterfly_ER_rcg");
elseif shape == "bird"
    bird_ER_rcg = ER;
    save('bird_er_rcg.mat', "bird_ER_rcg");
end    



%% construct plots
%*****************

% plot registered shape
%----------------------
f = figure('visible', 'off');
plot(X2(:,1), X2(:,2), Xreg(:,1), Xreg(:,2), 'LineStyle','none', 'Marker','.', 'MarkerSize',2);
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
title('Proposed RCG', 'FontSize',15);
% legend('target', 'registered', 'FontSize',10, 'Location','southwest');
grid on;
exportgraphics(f, strcat(shape, '_rcg.png'));



fprintf('\n--- Robust Riemannian CG Registration (SE(2)) ---\n');
fprintf('Iterations: %d\n', info.k);
fprintf('Final RMS Error: %.5f\n', info.ER(end));