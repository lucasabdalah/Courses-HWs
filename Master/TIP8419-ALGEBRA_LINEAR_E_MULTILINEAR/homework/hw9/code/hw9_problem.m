close all;
clearvars;


%% Part 2:
d.RMC = 1000;
d.N = 3; 
d.I = [2 3 4]; 
d.R  = 4;
d.SNR_dB = 0:5:30;
d.NMSE_dB = zeros(d.RMC, length(d.SNR_dB));

for rmc = 1:1:d.RMC
    for ii_snr = 1:1:length(d.SNR_dB)
        A0 = nd.randn_complex(d.I(1), d.R);
        B0 = nd.randn_complex(d.I(2), d.R);
        C0 = nd.randn_complex(d.I(3), d.R);
        X0 = nd.kr_(nd.kr_(A0, B0), C0);
        
        noise = addNoise(X0, d.SNR_dB(ii_snr), prod(d.I), d.R);
        X = X0 + noise;
        
        factors = nd.mlskrf(X, d.N, d.I);
        Xhat = nd.kr_(nd.kr_(factors{1},factors{2}),factors{3});
        [~, d.NMSE_dB(rmc, ii_snr)] = nd.nmse(X0, Xhat);
    end
    fprintf(' ----- (%2.0f, %2.0f) -------\n', rmc, ii_snr)
end

fprintf('---------------- \n')

d.meanNMSE = mean(d.NMSE_dB,1);
fprintf('Mean NMSE: %2.2f dB \n', d.meanNMSE);
% fprintf('Mean NMSE for %d MC rounds: %2.2f dB \n', d. RMC, d.meanNMSE);



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

% Save Data
% disp('Saving data...')
% save('hw9_problem_data.mat','-struct','d');
% disp('Saved sucessfully')


%% Local Function -------------------------------------------------------------
function Noise = addNoise(X0, SNRdB, M, N)
    beta = nd.randn_complex(M, N);
    alpha = frob(X0)/(10^(SNRdB/10)*frob(beta));
    Noise = alpha.*beta;
end