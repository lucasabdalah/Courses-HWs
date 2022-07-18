%% [TIP8419 - Algebra Linear e Multilinear] - Homework 1
% by Lucas Abdalah
% ----------------------------
% Problem
% hw1_problem3.m
% Author: Lucas Abdalah
%

%% General Setup
data.MC = 500;
data.N = 2.^(1:7);
data.time_method_author= zeros(data.MC, length(data.N));
A = cell(data.MC, length(data.N));
B = cell(data.MC, length(data.N));

tStart = tic;
%% Matrix Computation
for ii = 1:1:length(data.N)
    A{ii} = complex(randn(data.N(ii), data.N(ii)), randn(data.N(ii),data.N(ii)));
    B{ii} = complex(randn(data.N(ii), data.N(ii)), randn(data.N(ii),data.N(ii)));
end


for RMC = 1:1:data.MC
    fprintf('.\n.\n.\nMonte Carlo Step: %2.0f \n- - - - - - - - - - -\n', RMC)

    for ii = 1:1:length(data.N)
        fprintf('Matrix Dimension Step: %2.0f \n', data.N(ii))
        [data.time_method_author(RMC,ii)] = method_author(A{ii},B{ii});
    end

    fprintf('.\n.\n.\n Elapsed Time: %2.1fs \n', sum(... 
        data.time_method_author(:), 1))
end    


data.tEnd = toc(tStart); % pair tic toc
fprintf('.\n.\n.\nTotal Time: %2.1fs \n', data.tEnd)


disp('Saving data...')
save('hw1_problem3_data.mat','-struct','data');
disp('Saved sucessfully')


%% ----------------------------------------------------------------------------
function C = kr_prod(A, B)
    
    N = size(A,2);
    if N == size(B,2)
        P = reshape(A,1,[],N);
        Q = reshape(B,[],1,N);
        C = P.*Q;
        C = reshape(C,[],N);
    else
        error('number of columns should be equal')
    end
end

function [time, C] = method_author(A, B)
    tic;
    C = kr_prod(A, B);
    time = toc;
end
