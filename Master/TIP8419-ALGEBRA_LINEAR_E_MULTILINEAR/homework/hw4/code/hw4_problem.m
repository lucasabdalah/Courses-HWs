close all;
clearvars;


%% Part 2:
RMC = 1000;
d1.SNR_dB = 0:5:30;
d1.I = 2; d1.P = 3;
d1.J = 4; d1.Q = 5;
d1.NMSE_dB = zeros(RMC, length(d1.SNR_dB));

for rmc = 1:1:RMC
    for ii_snr = 1:1:length(d1.SNR_dB)
        A0 = nd.randn_complex(d1.I, d1.P);
        B0 = nd.randn_complex(d1.J, d1.Q);
        X0 = nd.kron_(A0, B0);
        
        noise = addNoise(X0, d1.SNR_dB(ii_snr), d1.I*d1.J, d1.P*d1.Q);
        X = X0 + noise;
        
        [Ahat,Bhat] = nd.lskronf(X, d1.I, d1.P, d1.J, d1.Q);
        Xhat = nd.kron_(Ahat, Bhat);
        [~, d1.NMSE_dB(rmc, ii_snr)] = nd.nmse(X0, Xhat);
    end
    fprintf(' ----- (%2.0f, %2.0f) -------\n', rmc, ii_snr)
end

fprintf('---------------- \n')

%% Part 3: 
d2.SNR_dB = 0:5:30;
d2.I = 4; d2.P = 3;
d2.J = 8; d2.Q = 5;
d2.NMSE_dB = zeros(RMC, length(d2.SNR_dB));

for rmc = 1:1:RMC
    for ii_snr = 1:1:length(d2.SNR_dB)
        A0 = nd.randn_complex(d2.I, d2.P);
        B0 = nd.randn_complex(d2.J, d2.Q);
        X0 = nd.kron_(A0, B0);
        
        noise = addNoise(X0, d2.SNR_dB(ii_snr), d2.I*d2.J, d2.P*d2.Q);
        X = X0 + noise;
    
        [Ahat,Bhat] = nd.lskronf(X, d2.I, d2.P, d2.J, d2.Q);
        Xhat = nd.kron_(Ahat, Bhat);
        [~, d2.NMSE_dB(rmc, ii_snr)] = nd.nmse(X0, Xhat);
    end
    fprintf(' ----- (%2.0f, %2.0f) -------\n', rmc, ii_snr)
end

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
axis tight


%% Save Data
disp('Saving data...')
save('hw4_problem_data1.mat','-struct','d1');
save('hw4_problem_data2.mat','-struct','d2');
disp('Saved sucessfully')


%% Local Function -------------------------------------------------------------
function Noise = addNoise(X0, SNRdB, M, N)
    beta = nd.randn_complex(M, N);
    alpha = frob(X0)/(10^(SNRdB/10)*frob(beta));
    Noise = alpha*beta;
end