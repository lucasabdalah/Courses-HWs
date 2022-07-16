% hw10_problem2.m

clc; pause(0.1)
% close all;
clearvars;
startTimeStamp = datestr(now,'HH:MM:SS.FFF');
%% Part 2:
d.RMC = 10;
d.N = 3; 
d.I = [5 5 5]; 
d.J = [5 5 5];
d.SNR_dB = 0:5:30;
d.NMSE_dB_HOSVD = zeros(d.RMC, length(d.SNR_dB));
d.NMSE_dB_HOOI = zeros(d.RMC, length(d.SNR_dB));

for rmc = 1:1:d.RMC
    for ii_snr = 1:1:length(d.SNR_dB)

        A0 = nd.randn_complex(d.I(1), d.J(1));
        B0 = nd.randn_complex(d.I(2), d.J(2));
        C0 = nd.randn_complex(d.I(3), d.J(3));
        X0 = nd.kron_(nd.kron_(A0, B0), C0);
    
        noise = addNoise(X0, d.SNR_dB(ii_snr), prod(d.I), prod(d.J));
        X = X0 + noise;
        % disp('oi')
        factors = nd.mlskronf(X, d.I, d.J, 'hosvd');
        Xhat = nd.kron_(nd.kron_(factors{1},factors{2}), factors{3});
        [~, d.NMSE_dB_HOSVD(rmc, ii_snr)] = nd.nmse(X0, Xhat);

        factors = nd.mlskronf(X, d.I, d.J, 'hooi');
        Xhat = nd.kron_(nd.kron_(factors{1},factors{2}),factors{3});
        [~, d.NMSE_dB_HOOI(rmc, ii_snr)] = nd.nmse(X0, Xhat);

    end
    dispIt(rmc, d.RMC, ii_snr);
end


fprintf('---------------- \n')

%% NMSE Output
d.meanNMSE_dB_HOSVD = mean(d.NMSE_dB_HOSVD,1);
d.meanNMSE_dB_HOOI = mean(d.NMSE_dB_HOOI,1);
fprintf('Mean NMSE for %d MC rounds:\n', d.RMC);
fprintf('HOSVD: %2.2f dB \n', d.meanNMSE_dB_HOSVD);
fprintf('HOOI: %2.2f dB \n', d.meanNMSE_dB_HOOI);
fprintf('---------------- \n')
fprintf('Mean Diff d1 vs d2: %2.2f dB \n', d.meanNMSE_dB_HOSVD - d.meanNMSE_dB_HOOI);
fprintf('Mean Diff: %2.2f dB \n', mean(d.meanNMSE_dB_HOSVD - d.meanNMSE_dB_HOOI, 2));

%% Figure Results
h_problem = figure;
plot(d.SNR_dB, d.meanNMSE_dB_HOSVD,...
    'Color', 'blue',...        
    'LineStyle', '-.',...
    'LineWidth', 1.0,...
    'Marker', 'o',...
    'MarkerFaceColor', 'blue',...
    'MarkerSize', 6);
hold on
plot(d.SNR_dB, d.meanNMSE_dB_HOOI,...
    'Color', 'red',...        
    'LineStyle', '--',...
    'LineWidth', 1.0,...
    'Marker', 's',...
    'MarkerFaceColor', 'red',...
    'MarkerSize', 5);
hold off
xticks(d.SNR_dB);
xlabel("SNR (dB)")
ylabel("NMSE (dB)")
legend(["HOSVD", "HOOI"], 'Location', 'Best')
legend boxoff
grid on
% axis tight

%% Save Data
disp('Saving data...')
fileName = 'hw10_problem_data2.mat';
% save(fileName,'-struct','d');
disp('Saved sucessfully');
endTimeStamp = datestr(now,'HH:MM:SS.FFF');
fprintf('Start message is sent at time %s\n', startTimeStamp);
fprintf('Finish message is sent at time %s\n', endTimeStamp);


%% Local Function -------------------------------------------------------------
function Noise = addNoise(X0, SNRdB, M, N)
    beta = nd.randn_complex(M, N);
    alpha = frob(X0)/(10^(SNRdB/10)*frob(beta));
    Noise = alpha.*beta;
end
function dispIt(mc, RMC, snr)
    fprintf(' ----- (%2.0f/%2.0f, %2.0f) -------\n', mc, RMC, snr);
end