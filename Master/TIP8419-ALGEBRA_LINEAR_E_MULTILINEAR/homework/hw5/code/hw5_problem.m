%% [TIP8419 - Algebra Linear e Multilinear] - Homework 5
% by Lucas Abdalah
% ----------------------------
% 
% hw5_problem.m
% Author: Lucas Abdalah
%
clearvars; close all;

%% Problem 2
RMC = 1000;
M = 3; N = 3; P = 3; Q = 3;
Xstruct = [M, N, P, Q];
R = 3;
d.Rvsrkp = 1:1:M*P;
d.NMSE_dB = zeros(RMC, M*P);

for rmc = 1:1:RMC
    for defRank = d.Rvsrkp
        X = modelkpsvd(R, Xstruct);
        [U,S,V,rkp] = nd.kpsvd(X, Xstruct);
        
        Xhat = buildkpsvd(U, S, V, defRank, Xstruct);
        [~, d.NMSE_dB(rmc, defRank)] = nd.nmse(X, Xhat);
        
        % fprintf(['\n.\n.\nOriginal Matrix vs KPSVD estimation\n' ...
        %         '-----------------------------------\n' ... 
        %         'rank\t|\tNMSE (dB)\t\n-----------------------------------\n']);
        %         fprintf('%2.0f\t|\t%2.2f\t\n', defRank, d.NMSE_dB(rmc, defRank));
    end
    fprintf(' ----- (%2.0f, %2.0f) -------\n', rmc, defRank)
end

fprintf('---------------- \n')

d.meanNMSE_dB = mean(d.NMSE_dB,1);
fprintf(' Mean NMSE: %2.2f dB \n', d.meanNMSE_dB);


%% Figure Results
h_problem = figure;
semilogy(d.Rvsrkp, d.meanNMSE_dB,...
    'Color', 'blue',...        
    'LineStyle', '--',...
    'LineWidth', 1.0,...
    'Marker', 'o',...
    'MarkerFaceColor', 'blue',...
    'MarkerSize', 5);
xticks(d.Rvsrkp);
xlabel("Rank")
ylabel("NMSE (dB)")
grid on
axis tight


%% Save Data
disp('Saving data...')
save('hw5_problem_data.mat','-struct','d');
disp('Saved sucessfully')

%% Local Function -------------------------------------------------------------
function X = modelkpsvd(R, Xstruct)
    % Using the outer product between vectors to control the rank of a matrix.
    X = complex(zeros(Xstruct(1)*Xstruct(3),Xstruct(2)*Xstruct(4)),0);
    for i = 1:R
        X = X + nd.randn_complex(Xstruct(1)*Xstruct(3), 1) * ...
            nd.randn_complex(1 ,Xstruct(2)*Xstruct(4));
    end
end


function Xhat = buildkpsvd(U, S, V, R, Xstruct)
    Xhat = complex(zeros(Xstruct(1)*Xstruct(3),Xstruct(2)*Xstruct(4)),0);
    Uk = cell(R);
    Vk = cell(R);

    for r = 1:R
        Uk{r} = reshape(U(:,r),[Xstruct(1) Xstruct(2)]);
        Vk{r}  = reshape(conj(V(:,r)), [Xstruct(3) Xstruct(4)]); 
        Xhat = Xhat + S(r,r)*nd.kron_(Uk{r},Vk{r});
    end
    Xhat = conj(Xhat);
end