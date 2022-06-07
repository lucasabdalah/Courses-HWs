%% [TIP8419 - Algebra Linear e Multilinear] - Homework 0 
% by Lucas Abdalah
%------------------------------------------------------------------------------
% script: hw0.m
% 2022/04/18 - v1

%% Clear Ambient
close all; % clc;


%% Problem 1 - item (a) ------------------------------------------------------------------
% hw0_problem1 % generate data
% disp('Saving data...')
% save('hw0_data.mat','-struct','data');
% disp('Saved sucessfully')
data = load('hw0_data.mat');


%% Plot results
N = data.N;
t1 = mean(data.time_a_method1, 1);
t2 = mean(data.time_a_method2, 1);

h = figure();
semilogy(N,t1,...
    'Color', 'blue',...        
    'LineStyle', '--',...
    'LineWidth', 1.0,...
    'Marker', 'o',...
    'MarkerFaceColor', 'blue',...
    'MarkerSize', 5);

hold on
semilogy(N,t2,...
    'Color', 'red',...        
    'LineStyle', '-',...
    'LineWidth', 1.0,...
    'Marker', 'square',...
    'MarkerFaceColor', 'red',...
    'MarkerSize', 5);

hold off
xticks(N);
xlabel('Matrix Dimension, N')
ylabel('Time (s)')
legend(["Method 1", "Method 2"], 'Location', 'Best') % legend(leg, 'Location', 'Northeastoutside')
legend boxoff
grid on
axis tight

savefig_tight(h, "hw0-problem1-a", "both")


%% Problem 1 - item (b) ------------------------------------------------------------------