%% [TIP8419 - Algebra Linear e Multilinear] - Homework 1
% by Lucas Abdalah
% ----------------------------
% Problem
% hw1.m
% Author: Lucas Abdalah
%

% clc;
clearvars; close all; 

%% hw1_problem1
problem1 = load('hw1_problem1_data.mat');

% Plot results
h1 = figure();
loglog(problem1.N, mean(problem1.time_method_mine, 1),...
    'Color', 'blue',...        
    'LineStyle', '--',...
    'LineWidth', 1.0,...
    'Marker', 'o',...
    'MarkerFaceColor', 'blue',...
    'MarkerSize', 5);

hold on
loglog(problem1.N,mean(problem1.time_method_matlab, 1),...
    'Color', 'red',...        
    'LineStyle', '-',...
    'LineWidth', 1.0,...
    'Marker', 'square',...
    'MarkerFaceColor', 'red',...
    'MarkerSize', 5);

hold off
xticks(problem1.N);
xlabel('Matrix Dimension, N')
ylabel('Time (s)')
legend(["Author", "Matlab"], 'Location', 'Best') % legend(leg, 'Location', 'Northeastoutside')
legend boxoff
grid on
axis tight

savefig_tight(h1, "figures/hw1-problem1", "both")

%% hw1_problem2
% data = load('hw1_problem2_data.mat');