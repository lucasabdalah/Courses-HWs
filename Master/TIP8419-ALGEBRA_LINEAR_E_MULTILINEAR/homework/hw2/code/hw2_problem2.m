%% [TIP8419 - Algebra Linear e Multilinear] - Homework 2
% Author: Lucas Abdalah
% ----------------------------
% Problem 2: Pseudo-inverse run time vs. rows
% hw2_problem2.m
%

clearvars; close all; 

%% General Setup
data.MC = 500;
data.N = 2:2:10;
data.time_1 = zeros(data.MC, length(data.N)); 

tStart = tic;

%% Monte Carlo Method
for RMC = 1:1:data.MC
    fprintf('.\n.\n.\nMonte Carlo Step: %2.0f \n- - - - - - - - - - -\n', RMC)

    for ii = 1:1:length(data.N)
            %% Compute each method time
            fprintf('Matrix Dimension Step: %2.0f \n', data.N(ii))
            data.time_1(RMC, ii) = method_1(data.N(ii));
    end

end    

data.tEnd = toc(tStart); % pair tic toc
fprintf('.\n.\n.\nTotal Time: %2.1fs \n', data.tEnd)

% disp('Saving data...')
% save('hw2_problem2_data.mat','-struct','data');
% disp('Saved sucessfully')

%% Local Functions ------------------------------------------------------------

function [time] = method_1(N)
    A = cell(N,1);
    I = 4; R = 2;

    for ii = 1:1:N
        A{ii} = complex(randn(I, R), randn(I,R));
    end 

    tic;
    for ii = 1:1:N
        if ii == 1
            C = A{ii};
        else
            C = nd.kr_(C, A{ii});
        end
    end
    time = toc;
end