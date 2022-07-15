%% [TIP8419 - Algebra Linear e Multilinear] - Homework 9
% by Lucas Abdalah
% ----------------------------
% 
% hw9.m
% Author: Lucas Abdalah
%

clc; pause(0.1)
clearvars;
close all;

%% Problem 1 ----------------------------
fprintf('------------------ Problem 1 ------------------------- \n')    

df = load('Practice_9_krf_matrix_3D.mat');

[I1, R] = size(df.A); [I2, ~] = size(df.B); [I3, ~] = size(df.C);
I_N = [I1, I2, I3];
[~, N_mode] = size(I_N);
X_size = size(df.X);

factors = nd.mlskrf(df.X, N_mode, I_N); 
Xhat = nd.kr_(nd.kr_(factors{1},factors{2}),factors{3});

[~, X_NMSE_dB] = nd.nmse(df.X, Xhat);
[~, A_NMSE_dB] = nd.nmse(df.A, factors{1});
[~, B_NMSE_dB] = nd.nmse(df.B, factors{2});
[~, C_NMSE_dB] = nd.nmse(df.C, factors{3});

fprintf('.\n.\nNMSE for\n')    
fprintf('---------------------------------------------- \n')    
fprintf('\tX and X_hat with MLSKRF: %2.2f dB \n', X_NMSE_dB);
fprintf('\tA and A_hat with MLSKRF: %2.2f dB \n', A_NMSE_dB);
fprintf('\tB and B_hat with MLSKRF: %2.2f dB \n', B_NMSE_dB);
fprintf('\tC and C_hat with MLSKRF: %2.2f dB \n', C_NMSE_dB);

fprintf('---------------------------------------------- \n');

X_scale = real(Xhat./df.X);
A_scale = factors{1}./df.A;
B_scale = factors{2}./df.B;
C_scale = factors{3}./df.C;

fprintf('.\n.\nScale factor for X and X_hat with MLSKRF\n');
fprintf('---------------------------------------------- \n')    
for jj = 1:1:X_size(2)
    for ii = 1:1:X_size(1)
    fprintf('\tX_hat./X(:, %g): %.2g \n', jj, X_scale(ii,jj));
    end
    % fprintf('\tX_hat./X(1:%g, %g): %.2g \n',X_size(1), jj, X_scale(ii,jj));
end

fprintf('.\n.\nScale factor for A and A_hat with MLSKRF\n');
fprintf('---------------------------------------------- \n')    
for jj = 1:1:R
    for ii = 1:1:I1
        fprintf('\tA_hat./A(:,%g): %.2g \n', jj, A_scale(ii,jj));
    end
    % fprintf('\tA_hat./A(1:%g, %g): %.2g \n', I1, jj, A_scale(ii,jj));
end

fprintf('.\n.\nScale factor for B and B_hat with MLSKRF\n');
fprintf('---------------------------------------------- \n')    
for jj = 1:1:R
    for ii = 1:1:I2
        fprintf('\tB_hat./B(:,%g): %.2g \n', jj, B_scale(ii,jj));
    end
    % fprintf('\tB_hat./B(1:%g, %g): %.2g \n', I2, jj, B_scale(ii,jj));
end

fprintf('.\n.\nScale factor for C and C_hat with MLSKRF\n');
fprintf('---------------------------------------------- \n')    
for jj = 1:1:R
    for ii = 1:1:I3
        fprintf('\tC_hat./C(:,%g): %.2g \n', jj, C_scale(ii,jj));
    end
    % fprintf('\tC_hat./C(1:%g, %g): %.2g \n', I3, jj, C_scale(ii,jj));
end

%% Part 2
fprintf('.\n.\nPart 2 ---------------------------------------------- \n')    

% hw9_problem

d = load('hw9_problem_data.mat');

fprintf(['SNR (dB) Mean NMSE (dB) for %d MC rounds\n' ...
    '-----------------------------------\n' ... 
    'SNR \t|\tNMSE \t\n----------------------------------\n'], d.RMC);
fprintf('%2.0f\t|\t%2.2f\t\n', [d.SNR_dB; d.meanNMSE]);


%% Figure Results
h_problem = figure;
plot(d.SNR_dB, d.meanNMSE,...
    'Color', 'blue',...        
    'LineStyle', '--',...
    'LineWidth', 1.0,...
    'Marker', 'o',...
    'MarkerFaceColor', 'blue',...
    'MarkerSize', 5);
xticks(d.SNR_dB);
xlabel("SNR (dB)")
ylabel("NMSE (dB)")
legend boxoff
grid on
% axis tight

%% Save Figures
% savefig_tight(h_problem, "figures/hw9-problem2", "both")