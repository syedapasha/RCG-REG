clear all; close all; clc;


% shape = "bunny";
shape = "dragon";

[X1, ~] = read_ply_xyz(strcat('../data/',shape,'-s.ply')); % source
[X2, ~] = read_ply_xyz(strcat('../data/',shape,'-t.ply')); % target


if size(X1,2) == 3, X1 = X1'; end
if size(X2,2) == 3, X2 = X2'; end


X1 = X1 - mean(X1, 2);
X2 = X2 - mean(X2, 2);


% N = min([size(X1,2), size(X2,2), 5000]);  
% idx1 = round(linspace(1, size(X1,2), N));
% idx2 = round(linspace(1, size(X2,2), N));
% X1 = X1(:,idx1);
% X2 = X2(:,idx2);

tic
[R, t, errors, n_iter, Xreg] = rsd_icp_3d( ...
    X1, X2, ...
    'max_iter', 1000, ...     
    'min_iter', 3, ...        
    'eta', 1, ...            
    'tol', 1e-5);             
toc

if shape == "bunny"
    bunny_ER_rsda = errors;
    save('bunny_ER_rsda.mat', "bunny_ER_rsda");
elseif shape == "dragon"
    dragon_ER_rsda = errors;
    save('dragon_ER_rsda.mat', "dragon_ER_rsda");
end


%% plot registered shape
%----------------------
f = figure;
pcshow(Xreg', 'BackgroundColor','white');
ax = gca; ax.FontSize = 10; % ax.PlotBoxAspectRatio = [1 1 1];
title("RSDA", 'FontSize',15);
exportgraphics(f, strcat(shape, '_rsda.png'));


% figure('Name','Riemannian SD ICP-3D Registration','NumberTitle','off');
% subplot(2,2,1);
% scatter3(X1(1,:), X1(2,:), X1(3,:), 8, 'r', 'filled');
% title('Source Shape'); axis equal; grid on; view(3);
% 
% subplot(2,2,2);
% scatter3(X2(1,:), X2(2,:), X2(3,:), 8, 'g', 'filled');
% title('Target Shape'); axis equal; grid on; view(3);
% 
% subplot(2,2,3);
% scatter3(X2(1,:), X2(2,:), X2(3,:), 8, [0 0.6 0], 'filled', 'MarkerFaceAlpha',0.5); hold on;
% scatter3(Xreg(1,:), Xreg(2,:), Xreg(3,:), 8, [0 0 1], 'filled', 'MarkerFaceAlpha',0.7);
% title('RSD-ICP (Registered vs Target)'); axis equal; grid on; view(3); legend('Target','Registered');
% 
% subplot(2,2,4);
% plot(0:n_iter, errors, 'k--o', 'LineWidth', 1.5);
% xlabel('Iteration'); ylabel('RMS Error');
% title('Convergence'); grid on;

fprintf('\n--- Riemannian Steepest Descent ICP 3D Registration ---\n');
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