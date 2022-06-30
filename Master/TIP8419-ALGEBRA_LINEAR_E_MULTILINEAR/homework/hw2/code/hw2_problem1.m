%% [TIP8419 - Algebra Linear e Multilinear] - Homework 2
% by Lucas Abdalah
% ----------------------------
% Problem
% hw2_problem1.m
% Author: Lucas Abdalah
%

%% General Setup
data.MC = 500;
data.N = 2.^(1:7);
data.time_method_author= zeros(data.MC, length(data.N));
data.time_method_matlab = zeros(data.MC, length(data.N));
A = cell(data.MC, length(data.N));
B = cell(data.MC, length(data.N)); 

tStart = tic;