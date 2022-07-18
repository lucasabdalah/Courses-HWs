%% [TIP8419 - Algebra Linear e Multilinear] - Seminar
% Author: Lucas Abdalah
% ----------------------------
% Partial reproduction of Figure 1: 
%   Alternating group lasso for block-term tensor decomposition with application 
%   to ECG source separation 
%   J. H. de M. Goulart et al., IEEE TSP, 2020
%
% plot_experiment1.m
%

% clear vars;
% clc; 
close all; 

filename = 'MC_scen1_rmc500.mat';

load(filename)

h1 = figure();

[fx, x] = ecdf(err_blocks_TL_SR(:));
[fy, y] = ecdf(err_blocks_TL_TR(:));

semilogx(x, fx, 'Color', 'blue', 'LineStyle', '-', 'LineWidth', 2.0);
hold on
semilogx(y, fy, 'Color', 'red', 'LineStyle', '-.', 'LineWidth', 2.0);
hold off
ylim([0 1])
xlim([2e-3 9.9e-1])
legend('BTD-GN (6, 6, 6)','BTD-GN (6, 5, 4)', 'Location', 'southeast')
legend boxoff
grid on
daspect([0.5 1 1])

savefig_tight(h1, "figures/seminar-figure1", "both")