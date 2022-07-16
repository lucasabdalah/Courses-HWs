%% [TIP8419 - Algebra Linear e Multilinear] - Homework 10
% by Lucas Abdalah
% ----------------------------
% 
% hw10.m
% Author: Lucas Abdalah
%

clc; pause(0.1)
clearvars;
close all;


%% Problem 1 ----------------------------
fprintf('------------------ Problem 1 ------------------------- \n')    

df = load('Practice_10_kronf_matrix_3D.mat');

[I1, J1] = size(df.A); [I2, J2] = size(df.B); [I3, J3] = size(df.C);
I = [I1 I2 I3]; J = [J1 J2 J3];

%% HOSVD version
fprintf('------------------ HOSVD ------------------------- \n');

factors = nd.mlskronf(df.X,I,J,'hosvd');
df.Xhat = nd.kron_(nd.kron_(factors{1},factors{2}),factors{3});

[~, X_NMSE_dB] = nd.nmse(df.X, df.Xhat);
[~, A_NMSE_dB] = nd.nmse(df.A, factors{1});
[~, B_NMSE_dB] = nd.nmse(df.B, factors{2});
[~, C_NMSE_dB] = nd.nmse(df.C, factors{3});

fprintf('.\n.\nNMSE for HOSVD\n');
fprintf('---------------------------------------------- \n');  
fprintf('\tX and X_hat with MLSKronF: %2.2f dB \n', X_NMSE_dB);
fprintf('\tA and A_hat with MLSKronF: %2.2f dB \n', A_NMSE_dB);
fprintf('\tB and B_hat with MLSKronF: %2.2f dB \n', B_NMSE_dB);
fprintf('\tC and C_hat with MLSKronF: %2.2f dB \n', C_NMSE_dB);

fprintf('---------------------------------------------- \n');

%% HOOI version
fprintf('------------------ HOOI ------------------------- \n');
factors = nd.mlskronf(df.X,I,J,'hooi');
df.Xhat = nd.kron_(nd.kron_(factors{1},factors{2}),factors{3});

[~, X_NMSE_dB] = nd.nmse(df.X, df.Xhat);
[~, A_NMSE_dB] = nd.nmse(df.A, factors{1});
[~, B_NMSE_dB] = nd.nmse(df.B, factors{2});
[~, C_NMSE_dB] = nd.nmse(df.C, factors{3});

fprintf('.\n.\nNMSE for HOOI\n');
fprintf('---------------------------------------------- \n');  
fprintf('\tX and X_hat with MLSKronF: %2.2f dB \n', X_NMSE_dB);
fprintf('\tA and A_hat with MLSKronF: %2.2f dB \n', A_NMSE_dB);
fprintf('\tB and B_hat with MLSKronF: %2.2f dB \n', B_NMSE_dB);
fprintf('\tC and C_hat with MLSKronF: %2.2f dB \n', C_NMSE_dB);

fprintf('---------------------------------------------- \n');

pause

d1 = load('hw10_problem_data1.mat');
d2 = load('hw10_problem_data2.mat');
d3 = load('hw10_problem_data3.mat');
d4 = load('hw10_problem_data4.mat');

%% Problem 1
h_problem1 = figure();
plot(d1.SNR_dB, d1.meanNMSE_dB_HOSVD,...
    'Color', 'blue',...        
    'LineStyle', '-.',...
    'LineWidth', 1.0,...
    'Marker', 'o',...
    'MarkerFaceColor', 'blue',...
    'MarkerSize', 6);
hold on
plot(d1.SNR_dB, d1.meanNMSE_dB_HOOI,...
    'Color', 'red',...        
    'LineStyle', '--',...
    'LineWidth', 1.0,...
    'Marker', 's',...
    'MarkerFaceColor', 'red',...
    'MarkerSize', 5);
hold off
xticks(d1.SNR_dB);
xlabel("SNR (dB)")
ylabel("NMSE (dB)")
legend(["HOSVD", "HOOI"], 'Location', 'northeast')
legend boxoff
grid on
axis tight
pbaspect([1 0.5 1])

savefig_tight(h_problem1, "figures/hw10-problem1", "both")


%% Problem 2
h_problem2 = figure();
plot(d2.SNR_dB, d2.meanNMSE_dB_HOSVD,...
    'Color', 'blue',...        
    'LineStyle', '-.',...
    'LineWidth', 1.0,...
    'Marker', 'o',...
    'MarkerFaceColor', 'blue',...
    'MarkerSize', 6);
hold on
plot(d2.SNR_dB, d2.meanNMSE_dB_HOOI,...
    'Color', 'red',...        
    'LineStyle', '--',...
    'LineWidth', 1.0,...
    'Marker', 's',...
    'MarkerFaceColor', 'red',...
    'MarkerSize', 5);
hold off
xticks(d2.SNR_dB);
xlabel("SNR (dB)")
ylabel("NMSE (dB)")
legend(["HOSVD", "HOOI"], 'Location', 'northwest')
legend boxoff
grid on
axis tight
pbaspect([1 0.5 1])

savefig_tight(h_problem2, "figures/hw10-problem2", "both")


%% Problem 3
h_problem3 = figure();
plot(d3.SNR_dB, d3.meanNMSE_dB_HOSVD,...
    'Color', 'blue',...        
    'LineStyle', '-.',...
    'LineWidth', 1.0,...
    'Marker', 'o',...
    'MarkerFaceColor', 'blue',...
    'MarkerSize', 6);
hold on
plot(d3.SNR_dB, d3.meanNMSE_dB_HOOI,...
    'Color', 'red',...        
    'LineStyle', '--',...
    'LineWidth', 1.0,...
    'Marker', 's',...
    'MarkerFaceColor', 'red',...
    'MarkerSize', 5);
hold off
xticks(d3.SNR_dB);
xlabel("SNR (dB)")
ylabel("NMSE (dB)")
legend(["HOSVD", "HOOI"], 'Location', 'northeast')
legend boxoff
grid on
axis tight
pbaspect([1 0.5 1])

savefig_tight(h_problem3, "figures/hw10-problem3", "both")


%% Problem 4
h_problem4 = figure();
plot(d4.SNR_dB, d4.meanNMSE_dB_HOSVD,...
    'Color', 'blue',...        
    'LineStyle', '-.',...
    'LineWidth', 1.0,...
    'Marker', 'o',...
    'MarkerFaceColor', 'blue',...
    'MarkerSize', 6);
hold on
plot(d4.SNR_dB, d4.meanNMSE_dB_HOOI,...
    'Color', 'red',...        
    'LineStyle', '--',...
    'LineWidth', 1.0,...
    'Marker', 's',...
    'MarkerFaceColor', 'red',...
    'MarkerSize', 5);
hold off
xticks(d4.SNR_dB);
xlabel("SNR (dB)")
ylabel("NMSE (dB)")
pbaspect([1 0.5 1])
legend(["HOSVD", "HOOI"], 'Location', 'northeast')
legend boxoff
grid on
axis tight

savefig_tight(h_problem4, "figures/hw10-problem4", "both")