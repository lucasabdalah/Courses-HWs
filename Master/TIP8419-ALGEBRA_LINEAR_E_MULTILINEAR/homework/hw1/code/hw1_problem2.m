%% [TIP8419 - Algebra Linear e Multilinear] - Homework 1
% by Lucas Abdalah
% ----------------------------
% Problem 2: Kronecker Product
% hw1_problem2.m
% Author: Lucas Abdalah
%

%% General Setup
data.MC = 500;
data.N = 2.^(1:7);
data.time_method_author= zeros(data.MC, length(data.N));
data.time_method_matlab = zeros(data.MC, length(data.N));
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
        [data.time_method_matlab(RMC,ii)] = method_matlab(A{ii},B{ii});
    end

    fprintf('.\n.\n.\n Elapsed Time: %2.1fs \n', sum([ ... 
        data.time_method_author(:);  data.time_method_matlab(:)], 1))
end    


data.tEnd = toc(tStart); % pair tic toc
fprintf('.\n.\n.\nTotal Time: %2.1fs \n', data.tEnd)


disp('Saving data...')
save('hw1_problem2_data.mat','-struct','data');
disp('Saved sucessfully')


%% ----------------------------------------------------------------------------

function C = kron_prod(A, B)
%kron_prod - Description
%
% Syntax: C = kron_prod(A, B)
%
% Long description

[M_rows, N_columns] = size(B);
C = repelem(A, M_rows, N_columns).*repmat(B, size(A));

end


% Author
function [time] = method_author(A, B)
    tic;
    kron_prod(A, B);
    time = toc;
end


% Matlab
function [time] = method_matlab(A, B)
    tic;
    kron(A,B);
    time = toc;
end