clear all; close all; 

shape = "lizard";
% shape = "dog";
% shape = "butterfly";
% shape = "bird";

%% Load source and target shapes
%*******************************
srcFile = strcat('../data/',shape,'-s.mat');
tgtFile = strcat('../data/',shape,'-t.mat');

srcMat = load(srcFile);     
tgtMat = load(tgtFile);     

if isfield(srcMat, 'X1')
    X1 = srcMat.X1;
elseif isfield(srcMat, 'source')
    X1 = srcMat.source;
elseif isfield(srcMat, 'contour')
    X1 = srcMat.contour;
else
    error(strcat('Source variable not found in ', srcFile));
end


if isfield(tgtMat, 'X2')
    X2 = tgtMat.X2;
elseif isfield(tgtMat, 'target')
    X2 = tgtMat.target;
elseif isfield(tgtMat, 'contour')
    X2 = tgtMat.contour;
else
    error(strcat('Target variable not found in ', srcFile));
end


% if size(X1,2) ~= 2
%     X1 = X1';
% end
% if size(X2,2) ~= 2
%     X2 = X2';
% end
if size(X1,1) ~= 2
    X1 = X1';
end
if size(X2,1) ~= 2
    X2 = X2';
end


% mean centering
X1 = X1 - mean(X1, 2);
X2 = X2 - mean(X2, 2);


N = min([size(X1,2), size(X2,2), 2000]);
idx1 = round(linspace(1, size(X1,2), N));
idx2 = round(linspace(1, size(X2,2), N));
X1 = X1(:,idx1);
X2 = X2(:,idx2);


[R, t, ER, n_iter, Xreg] = riemannian_sd_icp_2d( ...
    X1, X2, ...
    'max_iter', 5000, ...     
    'min_iter', 3, ...        
    'eta', 0.5, ...            
    'tol', 1e-5);             

if shape == "cattle"
    cattle_ER_rsda = ER;
    save('cattle_ER_rsda.mat', "cattle_ER_rsda");
elseif shape == "lizard"
    lizard_ER_rsda = ER;
    save('lizard_ER_rsda.mat', "lizard_ER_rsda");
elseif shape == "dog"
    dog_ER_rsda = ER;
    save('dog_ER_rsda.mat', "dog_ER_rsda");
end    



%% plot registered shape
%-----------------------
f = figure('visible', 'off');
plot(X2(1,:), X2(2,:), Xreg(1,:), Xreg(2,:), 'LineStyle','none',  'Marker','.', 'MarkerSize',2);
% xlim([-1,1]); ylim([-1,1]); 
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
title('RSDA', 'FontSize',15);
legend('target', 'registered', 'FontSize',10, 'Location','southwest');
grid on;
exportgraphics(f, strcat(shape, '_rsda.png'));


% figure('Name','Riemannian SD ICP-2D Registration','NumberTitle','off');
% subplot(2,2,1);
% scatter(X1(1,:), X1(2,:), 8, 'r', 'filled');
% title('Source Shape'); axis equal; grid on;
% 
% subplot(2,2,2);
% scatter(X2(1,:), X2(2,:), 8, 'g', 'filled');
% title('Target Shape'); axis equal; grid on;
% 
% subplot(2,2,3);
% hold on;
% scatter(X1(1,:), X1(2,:), 8, 'r', 'filled');         
% scatter(X2(1,:), X2(2,:), 8, 'g', 'filled');        
% scatter(Xreg(1,:), Xreg(2,:), 8, 'b', 'filled');    
% title('Source, Target, Registered Shape');
% axis equal; grid on; 
% legend('Source','Target','Registered');
% hold off;
% 
% subplot(2,2,4);
% plot(0:n_iter, errors, 'k--o', 'LineWidth', 1.5);
% xlabel('Iteration'); ylabel('RMS Error');
% title('Convergence'); grid on;


fprintf('\n--- Riemannian Steepest Descent ICP 2D Registration ---\n');
fprintf('Iterations: %d\n', n_iter);
fprintf('Initial RMS Error: %.6f\n', ER(1));
fprintf('Final RMS Error: %.6f\n', ER(end));