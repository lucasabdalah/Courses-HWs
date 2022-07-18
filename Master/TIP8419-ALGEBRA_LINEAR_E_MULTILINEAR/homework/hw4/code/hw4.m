%% [TIP8419 - Algebra Linear e Multilinear] - Homework 4
% by Lucas Abdalah
% ----------------------------
% 
% hw4.m
% Author: Lucas Abdalah
%

% clearvars; close all;
pauseTime = 5;

%% Part 1 - Experiment
fprintf('------------------ Experiment ---------------------------- \n')    

Asize = [4, 2];
Bsize = [6, 3];

A = nd.randn_complex(Asize(1), Asize(2));
B = nd.randn_complex(Bsize(1), Bsize(2));
X = nd.kron_(A,B);

[Ahat,Bhat] = nd.lskronf(X, Asize(1), Asize(2), Bsize(1), Bsize(2));
Xhat = nd.kron_(Ahat,Bhat);

[~, X_NMSE_dB] = nd.nmse(X, Xhat);
[~, A_NMSE_dB] = nd.nmse(A, Ahat);
[~, B_NMSE_dB] = nd.nmse(B, Bhat);

fprintf('.\n.\nNMSE for\n')    
fprintf('---------------------------------------------- \n')    
fprintf('\tX and X_hat with LSKRONF: %2.2f dB \n', X_NMSE_dB);
fprintf('\tA and A_hat with LSKRONF: %2.2f dB \n', A_NMSE_dB);
fprintf('\tB and B_hat with LSKRONF: %2.2f dB \n', B_NMSE_dB);
fprintf('---------------------------------------------- \n')    

A_scale = Ahat./A;
B_scale = Bhat./B;

fprintf('.\n.\nScale factor for A and A_hat with LSKRONF\n');
fprintf('---------------------------------------------- \n')    
for jj = 1:1:Asize(2)
    for ii = 1:1:Asize(1)
        fprintf('\tA_hat./A(:,%g): %.2g \n', jj, A_scale(ii,jj));
    end
end

fprintf('.\n.\nScale factor for B and B_hat with LSKRONF\n');
fprintf('---------------------------------------------- \n')    
for jj = 1:1:Bsize(2)
    for ii = 1:1:Bsize(1)
        fprintf('\tB_hat./B(:,%g): %.2g \n', jj, B_scale(ii,jj));
    end
end

pause(pauseTime); clc; pause(pauseTime);

fprintf('------------------- Validation --------------------------- \n')


%% Part 1 - Validation
load('Practice_4_kronf_matrix.mat');
[iX, jX] = size(X);
[iA, jA] = size(A);
[iB, jB] = size(B);

[Ahat,Bhat] = nd.lskronf(X, iA, jA, iB, jB);
Xhat = nd.kron_(Ahat,Bhat);

[~, X_NMSE_dB] = nd.nmse(X, Xhat);
[~, A_NMSE_dB] = nd.nmse(A, Ahat);
[~, B_NMSE_dB] = nd.nmse(B, Bhat);

fprintf('.\n.\nNMSE for\n')    
fprintf('---------------------------------------------- \n')    
fprintf('\tX and X_hat with LSKRONF: %2.2f dB \n', X_NMSE_dB);
fprintf('\tA and A_hat with LSKRONF: %2.2f dB \n', A_NMSE_dB);
fprintf('\tB and B_hat with LSKRONF: %2.2f dB \n', B_NMSE_dB);
fprintf('---------------------------------------------- \n')    

A_scale = Ahat./A;
B_scale = Bhat./B;

fprintf('.\n.\nScale factor for A and A_hat with LSKRONF\n');
fprintf('---------------------------------------------- \n')    
for jj = 1:1:jA
    for ii = 1:1:iA
        fprintf('\tA_hat./A(:,%g): %.2g \n', jj, A_scale(ii,jj));
    end
end

fprintf('.\n.\nScale factor for B and B_hat with LSKRONF\n');
fprintf('---------------------------------------------- \n')    
for jj = 1:1:jB
    for ii = 1:1:iB
        fprintf('\tB_hat./B(:,%g): %.2g \n', jj, B_scale(ii,jj));
    end
end


%% Part 2
fprintf('.\n.\nPart 2 ---------------------------------------------- \n')    

hw4_problem

d1 = load('hw4_problem_data1.mat');
d2 = load('hw4_problem_data2.mat');

d1.meanNMSE = mean(d1.NMSE_dB,1);
d2.meanNMSE = mean(d2.NMSE_dB,1);
fprintf(' Mean Diff d1 vs d2: %2.2f dB \n', d1.meanNMSE - d2.meanNMSE);
fprintf(' Mean Diff: %2.2f dB \n', mean(d1.meanNMSE - d2.meanNMSE, 2));


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
legend({sprintf('(I, P) (J, Q) = (%2.0f, %2.0f) (%2.0f, %2.0f)', ...
                d1.I, d1.P, d1.J, d1.Q), ... 
        sprintf('(I, P) (J, Q) = (%2.0f, %2.0f) (%2.0f, %2.0f)', ...
                d2.I, d2.P, d2.J, d2.Q)}, ...
        'Location', 'Best')
legend boxoff
grid on
% axis tight


%% Save Figures
% savefig_tight(h_problem, "figures/hw4-problem1", "both")