%% [TIP8419 - Algebra Linear e Multilinear] - Homework 0 
% by Lucas Abdalah
%------------------------------------------------------------------------------
% script: hw0.m
% 2022/04/18 - v1

%% Problem 1 ------------------------------------------------------------------

%% generate data

% item (a) 
% data = problem_1;
% disp('Saving data...')
% % save('hw0_problem1_data.mat', data)
% save('newstruct.mat','-struct','data');
% disp('Saved sucessfully')

%  Compute the mean for each method
% load('newstruct.mat')

%% plot results
close all
clc;
t1 = mean(time_method1, 1);
t2 = mean(time_method2, 1);
figure();
semilogy(N,t1);
hold on
semilogy(N,t2);
hold off
xticks(N);
xlabel('Matrix Dimension, N')
ylabel('Time (s)')
legend("Method 1", "Method 2")


% item (b) 

%% Problem 2 ------------------------------------------------------------------

