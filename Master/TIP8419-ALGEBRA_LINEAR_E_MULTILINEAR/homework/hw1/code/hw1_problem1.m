%% [TIP8419 - Algebra Linear e Multilinear] - Homework 1
% by Lucas Abdalah
% ----------------------------
% Problem
% hw1_problem1.m
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
save('hw1_problem1_data.mat','-struct','data');
disp('Saved sucessfully')

%% ----------------------------------------------------------------------------

function C = hadamard_prod(A, B)
%hadamard_prod - Description
%
% Syntax: C = hadamard_prod(A, B)
%
% Long description
if size(A) == size(B)
    % C = A.*B
    [Ai, Aj] = size(A);
    C = zeros(Ai, Aj);
    for ii = 1:1:Ai
        for jj = 1:1:Aj
            C(ii,jj) = A(ii, jj)*B(ii,jj);
        end
    end
else
    error('number of dimensions should be equal')
end

end


% Author
function [time, C] = method_author(A, B)
    tic;
    C = hadamard_prod(A, B);
    time = toc;
end


% Matlab
function [time, C] = method_matlab(A, B)
    tic;
    C = A.*B;
    time = toc;
end