% hw11_problem.m

clc; pause(0.1)
close all;
clearvars;


%% Part 2:
d.RMC = 1000;
d.I = 8; d.J = 4; d.K = 5; d.R = 3;
d.SNR_dB = 0:5:30;
d.NMSE_dB_X = zeros(d.RMC, length(d.SNR_dB));
d.NMSE_dB_A = zeros(d.RMC, length(d.SNR_dB));
d.NMSE_dB_B = zeros(d.RMC, length(d.SNR_dB));
d.NMSE_dB_C = zeros(d.RMC, length(d.SNR_dB));

for rmc = 1:d.RMC
    for ii_snr = 1:length(d.SNR_dB)

        A = nd.randn_complex(d.I, d.R);
        B = nd.randn_complex(d.J, d.R);
        C = nd.randn_complex(d.K, d.R);
        X0 = nd.fold(A*(nd.kr_(C,B).'),[d.I d.J d.K],1);
        noise = addNoise(X0, d.SNR_dB(ii_snr), d.I, d.J, d.K);

        X = X0 + noise;
        
        [ia,ib,ic] = size(X);
        [Ahat, Bhat, Chat, error, it] = nd.als(X, d.R, eps, 5e2);
        Xhat = nd.fold(Ahat*(nd.kr_(Chat,Bhat).'), [d.I d.J d.K], 1);
    
        [~, d.NMSE_dB_X(rmc, ii_snr)] = nd.nmse(X0, Xhat);
        [~, d.NMSE_dB_A(rmc, ii_snr)] = nd.nmse(A, Ahat);
        [~, d.NMSE_dB_B(rmc, ii_snr)] = nd.nmse(B, Bhat);
        [~, d.NMSE_dB_C(rmc, ii_snr)] = nd.nmse(C, Chat);
    end
    dispIt(rmc, ii_snr);
end


fprintf('----------------------------------- \n')


%% NMSE Output
d.meanNMSE_dB_X = mean(d.NMSE_dB_X,1);
d.meanNMSE_dB_A = mean(d.NMSE_dB_A,1);
d.meanNMSE_dB_B = mean(d.NMSE_dB_B,1);
d.meanNMSE_dB_C = mean(d.NMSE_dB_C,1);
fprintf('Mean NMSE for %d MC rounds:\n', d.RMC);


fprintf(['SNR (dB) Mean NMSE (dB) for %d MC rounds\n' ...
    '-----------------------------------\n' ... 
    'SNR \t|\tNMSE \t\n----------------------------------\n'], d.RMC);

fprintf('\t X \n-----------------------------------\n');
fprintf('%2.0f\t|\t%2.2f\t\n', [d.SNR_dB; d.meanNMSE_dB_X]);
fprintf('-----------------------------------\n')

fprintf('\t A \n-----------------------------------\n');
fprintf('%2.0f\t|\t%2.2f\t\n', [d.SNR_dB; d.meanNMSE_dB_A]);
fprintf('-----------------------------------\n')

fprintf('\t B \n-----------------------------------\n');;
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


%% Save Data
disp('Saving data...')
% save('hw11_problem_data.mat','-struct','d');
disp('Saved sucessfully')


%% Local Function -------------------------------------------------------------
function Noise = addNoise(X0, SNRdB, M, varargin)
    beta = nd.randn_complex(M, varargin{:});
    alpha = frob(X0)/(10^(SNRdB/10)*frob(beta));
    Noise = alpha.*beta;
end
function dispIt(mc, snr)
    fprintf(' ----- (%2.0f, %2.0f) -------\n', mc, snr);
end