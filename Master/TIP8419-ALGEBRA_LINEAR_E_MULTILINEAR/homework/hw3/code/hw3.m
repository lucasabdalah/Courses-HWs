%% [TIP8419 - Algebra Linear e Multilinear] - Homework 3
% by Lucas Abdalah
% ----------------------------
% 
% hw3.m
% Author: Lucas Abdalah
%

clearvars; close all;

%% Part 1 
load('Practice_3_krf_matrix.mat');
[iX, jX] = size(X);
[iA, jA] = size(A);
[iB, jB] = size(B);

[Ahat,Bhat] = nd.lskrf(X, iA, iB);
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
fprintf('.\n.\nPart 2 ---------------------------------------------- \n')    

hw3_problem

d1 = load('hw3_problem_data1.mat');
d2 = load('hw3_problem_data2.mat');

fprintf('\tMean Diff d1 vs d2: %2.2f dB \n', d1.meanNMSE - d2.meanNMSE);
fprintf('\tMean Diff: %2.2f dB \n', mean(d1.meanNMSE - d2.meanNMSE, 2));


%% Figure Results
h_problem = figure;
plot(d1.SNR_dB, d1.meanNMSE,...
    'Color', 'blue',...        
    'LineStyle', '--',...
    'LineWidth', 1.0,...
    'Marker', 'o',...
    'MarkerFaceColor', 'blue',...
    'MarkerSize', 5);
hold on
plot(d2.SNR_dB, d2.meanNMSE,...
    'Color', 'red',...        
    'LineStyle', '-',...
    'LineWidth', 1.0,...
    'Marker', 's',...
    'MarkerFaceColor', 'red',...
    'MarkerSize', 5);
hold off
xticks(d1.SNR_dB);
xlabel("SNR (dB)")
ylabel("NMSE (dB)")
legend(["I = " + num2str(d1.I), "I = " + num2str(d2.I)], 'Location', 'Best')
legend boxoff
grid on
% axis tight


%% Save Figures
savefig_tight(h_problem, "figures/hw3-problem1", "both")