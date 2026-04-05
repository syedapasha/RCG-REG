clear all; close all; clc; 

% shape = "bunny";
shape = "dragon";

[X1, ~] = read_ply_xyz(strcat('../data/',shape,'-s.ply')); % source
[X2, ~] = read_ply_xyz(strcat('../data/',shape,'-t.ply')); % target


if size(X1,2) == 3, X1 = X1'; end
if size(X2,2) == 3, X2 = X2'; end


% N = min([size(X1,2), size(X2,2), 1000]);
% idx1 = round(linspace(1, size(X1,2), N));
% idx2 = round(linspace(1, size(X2,2), N));
% X1 = X1(:,idx1);
% X2 = X2(:,idx2);

tic
[R, t, errors, n_iter, Xreg] = ehl_icp_3d(X1, X2, ...
    'max_iter', 1000, 'min_iter', 5, 'eta', .9, 'mu', 0.5, 'tol', 1e-6);
toc

if shape == "bunny"
    bunny_ER_ehl_icp = errors;
    save('bunny_ER_ehl_icp.mat', "bunny_ER_ehl_icp");
elseif shape == "dragon"
    dragon_ER_ehl_icp = errors;
    save('dragon_ER_ehl_icp.mat', "dragon_ER_ehl_icp");
end


all_points = [X1, X2, Xreg];
limx = [min(all_points(1,:)), max(all_points(1,:))];
limy = [min(all_points(2,:)), max(all_points(2,:))];
limz = [min(all_points(3,:)), max(all_points(3,:))];


%% plot registered shape
%----------------------
f = figure;
pcshow(Xreg', 'BackgroundColor','white');
ax = gca; ax.FontSize = 10; % ax.PlotBoxAspectRatio = [1 1 1];
title("EHL-ICP", 'FontSize',15);
exportgraphics(f, strcat(shape, '_ehl_icp.png'));


% figure('Name','EHL-ICP-3D Registration','NumberTitle','off');
% 
% 
% subplot(2,2,1);
% scatter3(X1(1,:), X1(2,:), X1(3,:), 8, 'r', 'filled');
% title('Source Shape');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% axis equal; grid on; view(3);
% xlim(limx); ylim(limy); zlim(limz);
% 
% 
% subplot(2,2,2);
% scatter3(X2(1,:), X2(2,:), X2(3,:), 8, 'g', 'filled');
% title('Target Shape');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% axis equal; grid on; view(3);
% xlim(limx); ylim(limy); zlim(limz);
% 
% % Registered vs Target
% subplot(2,2,3);
% scatter3(X2(1,:), X2(2,:), X2(3,:), 8, [0 0.6 0], 'filled', 'MarkerFaceAlpha',0.5); hold on;
% scatter3(Xreg(1,:), Xreg(2,:), Xreg(3,:), 8, [0 0 1], 'filled', 'MarkerFaceAlpha',0.7);
% title('Registered vs Target');
% xlabel('X'); ylabel('Y'); zlabel('Z');
% axis equal; grid on; view(3);
% legend('Target','Registered');
% xlim(limx); ylim(limy); zlim(limz);
% 
% % RMS Error Curve
% subplot(2,2,4);
% plot(0:n_iter, errors, 'k--o', 'LineWidth', 1.5);
% xlabel('Iteration');
% ylabel('RMS Error');
% title('Convergence (RMS Error)');
% grid on;

fprintf('\n--- EHL-ICP 3D Registration ---\n');
fprintf('Iterations: %d\n', n_iter);
fprintf('Final RMS Error: %.6f\n', errors(end));

% --- Helper function for reading PLY (same as before) ---
function [V, F] = read_ply_xyz(filename)
    fid = fopen(filename,'r');
    tline = '';
    n_vert = 0; n_face = 0;
    while ischar(tline)
        tline = fgetl(fid);
        if contains(tline, 'element vertex')
            n_vert = sscanf(tline, 'element vertex %d');
        elseif contains(tline, 'element face')
            n_face = sscanf(tline, 'element face %d');
        elseif strcmp(tline, 'end_header')
            break;
        end
    end
    V = fscanf(fid, '%f %f %f', [3 n_vert])';
    if n_face > 0
        F = fscanf(fid, '3 %d %d %d', [3 n_face])'+1;
    else
        F = [];
    end
    fclose(fid);
end