clear all; close all; clc; 
% --- Load Source and Target 3D Shapes from PLY files ---
% Assumes 'data_3D/mbunny.ply' and 'data_3D/hbunny.ply' exist in the current directory

% shape = "bunny";
% shape = "dragon";
shape = "armadillo";

% Read source and target point clouds
srcFile = fullfile('../data', strcat(shape,'-s.ply'));
tgtFile = fullfile('../data', strcat(shape,'-t.ply'));
srcPC = pcread(srcFile); % Requires Matlab's Computer Vision Toolbox
tgtPC = pcread(tgtFile);

X1 = srcPC.Location; % Nx3 source points
X2 = tgtPC.Location; % Mx3 target points


% --- Optional: Downsample for speed if too large ---
% maxN = 10000;
% if size(X1,1) > maxN
%     idx = randperm(size(X1,1), maxN);
%     X1 = X1(idx,:);
% end
% if size(X2,1) > maxN
%     idx = randperm(size(X2,1), maxN);
%     X2 = X2(idx,:);
% end


% --- Parameters ---
opts.max_iter = 500;
opts.xtol = 1e-8;
opts.gtol = 1e-8;
opts.ftol = 1e-8;
opts.alpha_max = 0.9;                       
opts.alpha_min = 1e-6;
opts.beta_type = 'DY';
opts.armijo_eta = 1e-2; %1e-3;
opts.armijo_gamma = 0.9;                   
opts.armijo_maxiter = 12;
opts.retraction_type = 'cayley'; % or 'exp'
opts.gamma = 0.85;

tic
[result, info] = cg_se3_registration(X1, X2, opts);
toc

Xreg = result.X1_reg;

if shape == "bunny"
    bunny_ER_rcg = info.ER;
    save('bunny_ER_rcg.mat', "bunny_ER_rcg");
elseif shape == "dragon"
    dragon_ER_rcg = info.ER;
    save('dragon_ER_rcg.mat', "dragon_ER_rcg");
else
    armadillo_ER_rcg = info.ER;
    save('armadillo_ER_rcg.mat', "armadillo_ER_rcg");
end



%% plot source and target
%------------------------
f = figure;
pcshow(X1, 'BackgroundColor','white');
ax = gca; ax.FontSize = 10; %ax.PlotBoxAspectRatio = [1 1 1];
title("Armadillo", 'FontSize',15);
exportgraphics(f, strcat(shape, '_src.png'));

f = figure;
pcshow(X2, 'BackgroundColor','white');
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
title("Armadillo Target", 'FontSize',15);
exportgraphics(f, strcat(shape, '_tgt.png'));


%% plot registered shape
%----------------------
f = figure;
pcshow(Xreg, 'BackgroundColor','white');
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
title("Proposed RCG", 'FontSize',15);
exportgraphics(f, strcat(shape, '_rcg.png'));



% % --- Plot Results ---

% 
% subplot(1,2,2);
% plot3(X2(:,1), X2(:,2), X2(:,3), 'g.', 'MarkerSize', 4); hold on;
% plot3(result.X1_reg(:,1), result.X1_reg(:,2), result.X1_reg(:,3), 'b.', 'MarkerSize', 4);
% title('Registered Source vs Target');
% axis equal; grid on; legend('Target','Registered Source');
% figure('Name','SE(3) Registration','NumberTitle','off');
% subplot(1,2,1);
% plot3(X1(:,1), X1(:,2), X1(:,3), 'r.', 'MarkerSize', 4);
% axis equal; grid on;
% subplot(1,2,2);
% plot3(X2(:,1), X2(:,2), X2(:,3), 'g.', 'MarkerSize', 4);
% axis equal; grid on;
% figure;
% plot(result.ER, 'k--o', 'LineWidth', 1.5);
% xlabel('Iteration');
% ylabel('RMS Error');
% title('Convergence Plot');
% grid on;

fprintf('\n--- Robust Nonmonotone Riemannian CG Registration (SE(3)) ---\n');
fprintf('Iterations: %d\n', info.k);
fprintf('Final RMS Error: %.5f\n', info.ER(end));