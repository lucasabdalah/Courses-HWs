%% [TIP8419 - Algebra Linear e Multilinear] - Homework 1
% by Lucas Abdalah
% ----------------------------
% Problem
% hw2.m
% Author: Lucas Abdalah
%

% clc;
clearvars; close all; 


%% hw2_problem1
problem1 = load('hw2_problem1_data.mat');

%% Plot results
h_problem1a = figure();
loglog(problem1.N, mean(problem1.time_1(:,:,1), 1),...
    'Color', 'blue',...        
    'LineStyle', '--',...
    'LineWidth', 1.0,...
    'Marker', 'o',...
    'MarkerFaceColor', 'blue',...
    'MarkerSize', 5);
hold on
loglog(problem1.N,mean(problem1.time_2(:,:,1), 1),...
    'Color', 'red',...        
    'LineStyle', '-',...
    'LineWidth', 1.0,...
    'Marker', 'square',...
    'MarkerFaceColor', 'red',...
    'MarkerSize', 5);
loglog(problem1.N,mean(problem1.time_3(:,:,1), 1),...
    'Color', 'green',...        
    'LineStyle', '-',...
    'LineWidth', 1.0,...
    'Marker', 'square',...
    'MarkerFaceColor', 'green',...
    'MarkerSize', 5);
hold off
xticks(problem1.N);
xlabel('Number of Rows, I (R=2)')
ylabel('Time (s)')
legend(["Method 1", "Method 2", "Method 3"], 'Location', 'Best') % legend(leg, 'Location', 'Northeastoutside')
legend boxoff
grid on
axis tight

savefig_tight(h_problem1a, "figures/hw2-problem1a", "both")

h_problem1b = figure();
loglog(problem1.N, mean(problem1.time_1(:,:,2), 1),...
    'Color', 'blue',...        
    'LineStyle', '--',...
    'LineWidth', 1.0,...
    'Marker', 'o',...
    'MarkerFaceColor', 'blue',...
    'MarkerSize', 5);
hold on
loglog(problem1.N,mean(problem1.time_2(:,:,2), 1),...
    'Color', 'red',...        
    'LineStyle', '-',...
    'LineWidth', 1.0,...
    'Marker', 'square',...
    'MarkerFaceColor', 'red',...
    'MarkerSize', 5);
loglog(problem1.N,mean(problem1.time_3(:,:,2), 1),...
    'Color', 'green',...        
    'LineStyle', '-',...
    'LineWidth', 1.0,...
    'Marker', 'square',...
    'MarkerFaceColor', 'green',...
    'MarkerSize', 5);
hold off
xticks(problem1.N);
xlabel('Number of Rows, I (R = 4)')
ylabel('Time (s)')
legend(["Method 1", "Method 2", "Method 3"], 'Location', 'Best') % legend(leg, 'Location', 'Northeastoutside')
legend boxoff
grid on
axis tight

savefig_tight(h_problem1b, "figures/hw2-problem1b", "both")

%% hw2_problem2
% problem2 = load('hw2_problem2_data.mat');
