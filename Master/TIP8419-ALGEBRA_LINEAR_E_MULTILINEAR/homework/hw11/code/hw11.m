%% [TIP8419 - Algebra Linear e Multilinear] - Homework 11
% by Lucas Abdalah
% ----------------------------
% 
% hw11.m
% Author: Lucas Abdalah
%

clc; pause(0.1)
clear;
close all;

%% Part 1 ------------------------------

df = load('Practice_9_cpd_tensor.mat');
[Ia,Ib,Ic] = size(df.tenX);
R = 3;

[Ahat, Bhat, Chat, error, it] = nd.als(df.tenX, R);

Xhat = nd.fold(Ahat*(nd.kr_(Chat,Bhat).'),[Ia Ib Ic],1);


%% NMSE Output
[~, X_NMSE_dB] = nd.nmse(df.tenX, Xhat);
[~, A_NMSE_dB] = nd.nmse(df.A, Ahat);
[~, B_NMSE_dB] = nd.nmse(df.B, Bhat);
[~, C_NMSE_dB] = nd.nmse(df.C, Chat);

fprintf('.\n.\n\n--------------------- Part 1 --------------------- \n\n');
fprintf('.\n.\nNMSE for ALS Validation\n');
fprintf('---------------------------------------------- \n');  
fprintf('\tX and X_hat: %2.2f dB \n', X_NMSE_dB);
fprintf('\tA and A_hat: %2.2f dB \n', A_NMSE_dB);
fprintf('\tB and B_hat: %2.2f dB \n', B_NMSE_dB);
fprintf('\tC and C_hat: %2.2f dB \n', C_NMSE_dB);

fprintf('---------------------------------------------- \n');


%% Figure Results
h_problem = figure;
plot(1:it, db(error),...
    'Color', 'blue',...        
    'LineStyle', '--',...
    'LineWidth', 1.0,...
    'Marker', 'o',...
    'MarkerFaceColor', 'blue',...
    'MarkerSize', 5);
xlabel("Iteration")
ylabel("Error (dB)")
legend('ALS$(\mathcal{X}, R)$','location','best', 'Interpreter','latex');
legend boxoff
grid on


%% Save Figures
% savefig_tight(h_problem, "figures/hw11-problem1", "both")


%% Part 2 ------------------------------
% hw11_problem.m
d = load('hw11_problem_data.mat');

fprintf('.\n.\n\n--------------------- Part 2 --------------------- \n\n');


%% NMSE Output
fprintf(['SNR (dB) Mean NMSE (dB) for %d MC rounds\n' ...
    '--------------------------------------------------\n\n' ... 
    'SNR \t|\tNMSE \t\n----------------------------------\n'], d.RMC);
fprintf('\t X \n-----------------------------------\n');
fprintf('%2.0f\t|\t%2.2f\t\n', [d.SNR_dB; d.meanNMSE_dB_X]);
fprintf('-----------------------------------\n')
fprintf('\t A \n-----------------------------------\n');
fprintf('%2.0f\t|\t%2.2f\t\n', [d.SNR_dB; d.meanNMSE_dB_A]);
fprintf('-----------------------------------\n')
fprintf('\t B \n-----------------------------------\n');
fprintf('%2.0f\t|\t%2.2f\t\n', [d.SNR_dB; d.meanNMSE_dB_B]);
fprintf('-----------------------------------\n')
fprintf('\t C \n-----------------------------------\n');
fprintf('%2.0f\t|\t%2.2f\t\n', [d.SNR_dB; d.meanNMSE_dB_C]);
fprintf('-----------------------------------\n')


%% Figure Results
h_problem1 = figure;
plot(d.SNR_dB, d.meanNMSE_dB_X,...
    'Color', 'black',...        
    'LineStyle', '-',...
    'LineWidth', 1.0,...
    'Marker', 'o',...
    'MarkerFaceColor', 'black',...
    'MarkerSize', 6);
hold on
plot(d.SNR_dB, d.meanNMSE_dB_A,...
    'Color', 'blue',...        
    'LineStyle', '--',...
    'LineWidth', 1.0,...
    'Marker', 's',...
    'MarkerFaceColor', 'blue',...
    'MarkerSize', 5);
plot(d.SNR_dB, d.meanNMSE_dB_B,...
    'Color', 'red',...        
    'LineStyle', '--',...
    'LineWidth', 1.0,...
    'Marker', 's',...
    'MarkerFaceColor', 'red',...
    'MarkerSize', 5);
plot(d.SNR_dB, d.meanNMSE_dB_C,...
    'Color', 'green',...        
    'LineStyle', '--',...
    'LineWidth', 1.0,...
    'Marker', 's',...
    'MarkerFaceColor', 'green',...
    'MarkerSize', 5);
hold off
xticks(d.SNR_dB);
xlabel("SNR (dB)")
ylabel("NMSE (dB)")
legend({'$\mathcal{{X}}$', '$\mathbf{{A}}$', '$\mathbf{{B}}$', '$\mathbf{{C}}$'}, 'location', 'best', 'Interpreter', 'latex');
legend boxoff
grid on


%% Save Figures
% savefig_tight(h_problem, "figures/hw11-problem2", "both")