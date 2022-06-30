%% [TIP8419 - Algebra Linear e Multilinear] - Homework 2
% Author: Lucas Abdalah
% ----------------------------
% Problem 1: Pseudo-inverse run time vs. rows
% hw2_problem1.m
%

clearvars; close all; 

%% General Setup
data.MC = 500;
data.N = 2.^(1:8);
data.R = [2, 4];
data.time_1 = zeros(data.MC, length(data.N), length(data.R));
data.time_2 = zeros(data.MC, length(data.N), length(data.R));
data.time_3 = zeros(data.MC, length(data.N), length(data.R));
A = cell(length(data.N), length(data.R));
B = cell(length(data.N), length(data.R)); 

tStart = tic;

%% Monte Carlo Method
for RMC = 1:1:data.MC
    fprintf('.\n.\n.\nMonte Carlo Step: %2.0f \n- - - - - - - - - - -\n', RMC)

    for ii = 1:1:length(data.N)
        I = data.N(ii);
        for jj = 1:1:length(data.R)
            R = data.R(jj);
            A{ii, jj} = complex(randn(I, R), randn(I, R));
            B{ii, jj} = randn(I, R);
            
            %% Compute each method time
            fprintf('Matrix Dimension Step: %2.0f \n', data.N(ii))
            data.time_1(RMC, ii, jj) = method_1(A{ii, jj},B{ii, jj});
            data.time_2(RMC, ii, jj) = method_2(A{ii, jj},B{ii, jj});
            data.time_3(RMC, ii, jj) = method_3(A{ii, jj},B{ii, jj});
        end
    end

end    

data.tEnd = toc(tStart); % pair tic toc
fprintf('.\n.\n.\nTotal Time: %2.1fs \n', data.tEnd)

disp('Saving data...')
save('hw2_problem1_data.mat','-struct','data');
disp('Saved sucessfully')

%% Local Functions ------------------------------------------------------------

function [time, C] = method_1(A, B)
    tic;
    X = nd.kr_(A,B);
    C = pinv(X);
    time = toc;
end


function [time, C] = method_2(A, B)
    tic;
    X = nd.kr_(A,B);
    C = inv(X.'*X)*X.';
    time = toc;
end


function [time, C] = method_3(A, B)
    tic;
    X = nd.kr_(A,B);
    C = inv(nd.hadamard_(A.'*A, B.'*B))*X.';
    time = toc;
end
