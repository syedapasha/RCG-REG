clear all; close all; clc;

% shape = "cattle";
% shape = "lizard";
shape = "dog";

% %% show convergence for different retractions
% load(strcat('rcg/',shape,'_ER_rcg_exp.mat'));
% ER_exp = ER;
% load(strcat('rcg/',shape,'_ER_rcg_Id.mat'));
% ER_Id = ER;
% load(strcat('rcg/',shape,'_ER_rcg_cay.mat'));
% ER_cay = ER;
% 
% f = figure('visible', 'off');
% 
% semilogy(ER_cay, 'Color','#0072BD', 'LineStyle','-', 'Marker','.'); hold on;
% semilogy(ER_exp, 'Color','#D95319', 'LineStyle','-', 'Marker','.');
% semilogy(ER_Id, 'Color','#EDB120', 'LineStyle','-', 'Marker','.');
% ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
% title('RMS Error', 'FontSize',15);
% legend('Cayley', 'Exponential', 'Identity', 'FontSize',10, 'Location','northeast');
% grid on;
% exportgraphics(f, 'retract.png');



%% compare convergence for different methods

load(strcat('icp/',shape,'_er_icp.mat'));
load(strcat('ehl_icp/',shape,'_er_ehl_icp.mat'));
load(strcat('rsda/',shape,'_er_rsda.mat'));
load(strcat('rcg/',shape,'_er_rcg.mat'));

% dog_ER_icp = dog_ER_icp(2:length(dog_ER_icp));


%% plot convergence behavior
f = figure('visible', 'off');

if shape == "cattle"
    plot(cattle_ER_icp, 'Color','#0072BD', 'LineStyle','-', 'Marker','.'); hold on;
    plot(cattle_ER_ehl_icp, 'Color','#D95319', 'LineStyle','-', 'Marker','.');
    plot(cattle_ER_rsda, 'Color','#EDB120', 'LineStyle','-', 'Marker','.');
    plot(cattle_ER_rcg, 'Color','#7E2F8E', 'LineStyle','-', 'Marker','.');
elseif shape == "lizard"
    plot(lizard_ER_icp, 'Color','#0072BD', 'LineStyle','-', 'Marker','.'); hold on;
    plot(lizard_ER_ehl_icp, 'Color','#D95319', 'LineStyle','-', 'Marker','.');
    plot(lizard_ER_rsda, 'Color','#EDB120', 'LineStyle','-', 'Marker','.');
    plot(lizard_ER_rcg, 'Color','#7E2F8E', 'LineStyle','-', 'Marker','.');
elseif shape == "dog"
    semilogy(dog_ER_icp, 'Color','#0072BD', 'LineStyle','-', 'Marker','.'); hold on;
    semilogy(dog_ER_ehl_icp, 'Color','#D95319', 'LineStyle','-', 'Marker','.');
    semilogy(dog_ER_rsda, 'Color','#EDB120', 'LineStyle','-', 'Marker','.');
    semilogy(dog_ER_rcg, 'Color','#7E2F8E', 'LineStyle','-', 'Marker','.');
end    
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
title('RMS Error', 'FontSize',15);
legend('ICP', 'EHL-ICP', 'RSDA', 'Proposed RCG', 'FontSize',10, 'Location','northeast');
grid on;
exportgraphics(f, strcat(shape,'_convergence.png'));
