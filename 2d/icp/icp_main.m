clear all; close all;

% shape = "lizard";
% shape = "dog";
% shape = "butterfly";
shape = "bird";


%% Load source and target shapes
%*******************************
src_data = load(strcat('../data/',shape,'-s.mat'));
tgt_data = load(strcat('../data/',shape,'-t.mat'));

% Find data points for each
src_names = fieldnames(src_data);
tgt_names = fieldnames(tgt_data);

src_pts = src_data.(src_names{1});
tgt_pts = tgt_data.(tgt_names{1});

% Ensure data are 2 x N
if size(src_pts,2) == 2
    src_pts = src_pts';
end
if size(tgt_pts,2) == 2
    tgt_pts = tgt_pts';
end

% Convert to 3 x N for ICP (z=0 for 2D)
src_pts3 = [ src_pts; zeros(1, size(src_pts,2)) ];
tgt_pts3 = [ tgt_pts; zeros(1, size(tgt_pts,2)) ];



%% ICP registration
%******************
tic
[TR, TT, ER, t] = icp(tgt_pts3, src_pts3, 30, 'Matching', 'bruteForce');
toc

% Transform source points by ICP result
src_icp = TR * src_pts3 + repmat(TT, 1, size(src_pts3,2));


if shape == "lizard"
    lizard_ER_icp = ER;
    save('lizard_er_icp.mat', "lizard_ER_icp");
elseif shape == "dog"
    dog_ER_icp = ER;
    save('dog_er_icp.mat', "dog_ER_icp");
elseif shape == "butterfly"
    butterfly_ER_icp = ER;
    save('butterfly_er_icp.mat', "butterfly_ER_icp");
elseif shape == "bird"
    bird_ER_icp = ER;
    save('bird_er_icp.mat', "bird_ER_icp");
end    




%% construct plots
%*****************

% plot source and target
f = figure('visible', 'off');
plot(src_pts3(1,:), src_pts3(2,:), tgt_pts3(1,:), tgt_pts3(2,:), 'LineStyle','none', 'Marker','.', 'MarkerSize',2);
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
title(strcat(shape," Test Set"), 'FontSize',15);
% legend('source', 'target', 'FontSize',10, 'Location','southwest');
grid on;
exportgraphics(f, strcat(shape, '_data.png'));


% plot registered shape
f = figure('visible', 'off');
plot(tgt_pts3(1,:), tgt_pts3(2,:), src_icp(1,:), src_icp(2,:),  'LineStyle','none',  'Marker','.', 'MarkerSize',2);
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
title('ICP', 'FontSize',15);
% legend('target', 'registered', 'FontSize',10, 'Location','southwest');
grid on;
exportgraphics(f, strcat(shape, '_icp.png'));
