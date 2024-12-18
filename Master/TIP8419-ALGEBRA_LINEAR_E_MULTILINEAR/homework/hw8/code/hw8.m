%% [TIP8419 - Algebra Linear e Multilinear] - Homework 8
% by Lucas Abdalah
% ----------------------------
% 
% hw8.m
% Author: Lucas Abdalah
%

clc; pause(0.1)
clear;
close all;

%% Problem 1 ----------------------------
df = load('hooi_test.mat');
maxIt = 100;
[df.tenS_hat,U_hat, it] = nd.hooi(df.tenX, 'full', maxIt);
df.tenX_hat = nd.N_mode(df.tenS_hat,U_hat);

f_ort = nd.sliceort(df.tenS_hat);
fprintf('\n--------------- The tensor slices are orthogonal if = zero -------------- \n'); 
fprintf('\tsum(Dot product) = %2.0f \n', f_ort);
fprintf('\n\n')

fprintf('\n--------------- Number of Iterations -------------- \n'); 
fprintf('\tit = %d \n', it);
fprintf('\n\n')

fprintf('\n--------------- NMSE between a given tensor X and estimation -------------- \n');
[~, tenX_NMSEdB] = nd.nmse(df.tenX, df.tenX_hat);
fprintf('\tNMSE: %2.2f dB\n', tenX_NMSEdB);
fprintf('\n\n')

fprintf('\n--------------- NMSE between a given tensor core S and its estimation ----- \n');
[~, tenS_NMSEdB] = nd.nmse(df.tenS, df.tenS_hat);
fprintf('\tNMSE: %2.2f dB\n', tenS_NMSEdB);
fprintf('\n\n')

fprintf('\n--------------- NMSE between the factor matrices U and their estimation ----- \n');
[~, U_NMSEdB(1)] = nd.nmse(df.U1, U_hat{1});
[~, U_NMSEdB(2)] = nd.nmse(df.U2, U_hat{2});
[~, U_NMSEdB(3)] = nd.nmse(df.U3, U_hat{3});
fprintf('\tNMSE between U%d and its estimation: %2.2f dB\n', [1:3; U_NMSEdB]);
fprintf('\n\n')


%% Problem 2 ----------------------------
I = [8 4 10]; J = [5 5 5];

X = nd.randn_complex(I);
Y = nd.randn_complex(J);

[S1,U1] = nd.hooi(X, 'trunc', maxIt);
[S2,U2] = nd.hooi(Y, 'trunc', maxIt);

R = size(S1);
P = size(S2);

Xhat = nd.N_mode(S1,U1,3);
Yhat = nd.N_mode(S2,U2,3);


fprintf('\n--------------- NMSE between a given tensor X and its estimation ----- \n');
[~, X_NMSEdB] = nd.nmse(X, Xhat);
fprintf('\tNMSE: %2.2f dB\n', X_NMSEdB);
fprintf('\n\n')

fprintf('\n--------------- NMSE between a given tensor Y and its estimation ----- \n');
[~, Y_NMSEdB] = nd.nmse(Y, Yhat);
fprintf('\tNMSE: %2.2f dB\n', Y_NMSEdB);
fprintf('\n\n')

fprintf('Tensor X multilinear rank: [%2.0f %2.0f %2.0f] \n', R)
fprintf('Tensor Y multilinear rank: [%2.0f %2.0f %2.0f] \n', P)