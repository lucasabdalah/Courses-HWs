%% [TIP8419 - Algebra Linear e Multilinear] - Homework 3
% by Lucas Abdalah
% ----------------------------
% 
% hw3.m
% Author: Lucas Abdalah
%


clearvars; close all;
% clc;

%% Part 1 

load('Practice_3_krf_matrix.mat');
[iX, jX] = size(X);
[iA, jA] = size(A);
[iB, jB] = size(B);

init.A0 = complex(zeros(iA,jA));
init.B0 = complex(zeros(iB,jB));

[Ahat,Bhat] = nd.lskrf(X, init);
Xhat = nd.kr_(Ahat,Bhat);

[X_NMSE, X_NMSE_dB] = nd.nmse(X, Xhat);
[A_NMSE, A_NMSE_dB] = nd.nmse(A, Ahat);
[B_NMSE, B_NMSE_dB] = nd.nmse(B, Bhat);

fprintf('.\n.\nNMSE for\n')    
fprintf('---------------------------------------------- \n')    
fprintf('\tX and X_hat with LSKRF: %2.2f dB \n', X_NMSE_dB);
fprintf('\tA and A_hat with LSKRF: %2.2f dB \n', A_NMSE_dB);
fprintf('\tB and B_hat with LSKRF: %2.2f dB \n', B_NMSE_dB);
fprintf('---------------------------------------------- \n')    

A_scale = Ahat./A;
B_scale = Bhat./B;

fprintf('.\n.\nScale factor for A and A_hat with LSKRF\n');
fprintf('---------------------------------------------- \n')    
for jj = 1:1:jA
    for ii = 1:1:iA
        fprintf('\tA_hat./A(:,%g): %.2g \n', jj, A_scale(ii,jj));
    end
end

fprintf('.\n.\nScale factor for B and B_hat with LSKRF\n');
fprintf('---------------------------------------------- \n')    
for jj = 1:1:jB
    for ii = 1:1:iB
        fprintf('\tB_hat./B(:,%g): %.2g \n', jj, B_scale(ii,jj));
    end
end


%% Part 2
pause
% RMC = 10;
% SNR_dB = 0:5:30;
% I = 10; J = 10; R = 4;


% X0 = nd.kr_(A, B);

% X = X0 + alpha*V;

% A = nd.randn_complex(I, J);
% B = nd.randn_complex(J, J);

% Xhat = tensor.mtx_prod_kr(Ahat,Bhat);

%% Local Function -------------------------------------------------------------