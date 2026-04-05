clear all; close all; clc;

% shape = "bunny";
% shape = "dragon";
shape = "armadillo";

%% show convergence for different retractions
% load('rcg/armadillo_ER_exp.mat');
% ER_exp = armadillo_ER_rcg;
% load('rcg/armadillo_ER_rcg.mat');
% ER_cay = armadillo_ER_rcg;
% 
% f = figure('visible', 'off');
% 
% semilogy(ER_cay, 'Color','#0072BD', 'LineStyle','-', 'Marker','.'); hold on;
% semilogy(ER_exp, 'Color','#D95319', 'LineStyle','-', 'Marker','.');
% ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
% title('RMS Error', 'FontSize',15);
% legend('Cayley', 'Exponential', 'FontSize',10, 'Location','northeast');
% grid on;
% exportgraphics(f, strcat(shape,'_retract.png'));



%% compare convergence for different methods

load(strcat('ehl_icp/',shape,'_ER_ehl_icp.mat'));
load(strcat('rsda/',shape,'_ER_rsda.mat'));
load(strcat('rcg/',shape,'_ER_rcg.mat'));


%% plot convergence behavior
f = figure('visible', 'off');

if shape == "bunny"
    semilogy(bunny_ER_ehl_icp, 'Color','#D95319', 'LineStyle','-', 'Marker','.'); hold on;
    semilogy(bunny_ER_rsda, 'Color','#EDB120', 'LineStyle','-', 'Marker','.');
    semilogy(bunny_ER_rcg, 'Color','#7E2F8E', 'LineStyle','-', 'Marker','.');
elseif shape == "dragon"
    semilogy(dragon_ER_ehl_icp, 'Color','#D95319', 'LineStyle','-', 'Marker','.'); hold on;
    semilogy(dragon_ER_rsda, 'Color','#EDB120', 'LineStyle','-', 'Marker','.');
    semilogy(dragon_ER_rcg, 'Color','#7E2F8E', 'LineStyle','-', 'Marker','.');
else
    semilogy(armadillo_ER_ehl_icp, 'Color','#D95319', 'LineStyle','-', 'Marker','.');  hold on;
    semilogy(armadillon_ER_rsda, 'Color','#EDB120', 'LineStyle','-', 'Marker','.');
    semilogy(armadillo_ER_rcg, 'Color','#7E2F8E', 'LineStyle','-', 'Marker','.');
end    
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
title('RMS Error', 'FontSize',15);
legend('EHL-ICP', 'RSDA', 'Proposed RCG', 'FontSize',10, 'Location','northeast');
grid on;
exportgraphics(f, strcat(shape,'_convergence.png'));
