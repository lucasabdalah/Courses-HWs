%% [TIP8419 - Algebra Linear e Multilinear] - Homework 5
% by Lucas Abdalah
% ----------------------------
% 
% hw5.m
% Author: Lucas Abdalah
%

clearvars; close all;

%% Part 1 -------------------------------
M = 3; N = 3; P = 3; Q = 3;
Xstruct = [M, N, P, Q];
R = 3;

fprintf('\n.\n.\n ---------------- Problem 1 ----------------\n.\n.\n')
X = modelkpsvd(R, Xstruct);
[U,S,V,rkp] = nd.kpsvd(X, Xstruct);

Xhat = buildkpsvd(U, S, V, rkp, Xstruct);
[~, ndNMSEdB] = nd.nmse(X, Xhat);
fprintf(['Original Matrix vs KPSVD estimation (full rank):\n' ...
        '\tNMSE = %2.2f dB \n'], ndNMSEdB);


%% Part 2 -----------------------------------
% hw5_problem

fprintf('\n.\n.\n ---------------- Problem 2 ----------------\n.\n.\n')
d = load('hw5_problem_data.mat');

fprintf(['Original Matrix vs KPSVD estimation\n' ...
                '-----------------------------------\n' ... 
                'rank\t|\tNMSE (dB)\t\n----------------------------------\n']);

fprintf('%2.0f\t|\t%2.2f\t\n', [d.Rvsrkp; d.meanNMSE_dB]);


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

% savefig_tight(h_problem, "figures/hw5-problem2", "both")


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