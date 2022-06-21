%% [TIP8419 - Algebra Linear e Multilinear] - Homework 0 
% by Lucas Abdalah
% ----------------------------
% Problema 1 - Item a: Kronecker Product and Matrix Inversion Cost
% hw0_problem1.m
% Author: Lucas Abdalah
%

log_write("start")

%% item (a)

%% General Setup
data1.MC = 100;
data1.N = 2.^(1:6);
data1.time_a_method1 = zeros(data1.MC, length(data1.N));
data1.time_a_method2 = zeros(data1.MC, length(data1.N));
A = cell(data1.MC, length(data1.N));
B = cell(data1.MC, length(data1.N)); 

tStart = tic;

%% Matrix Computation
for ii = 1:1:length(data1.N)
    A{ii} = complex(randn(data1.N(ii), data1.N(ii)), randn(data1.N(ii),data1.N(ii)));
    B{ii} = complex(randn(data1.N(ii), data1.N(ii)), randn(data1.N(ii),data1.N(ii)));
end


for RMC = 1:1:data1.MC
    fprintf('.\n.\n.\nMonte Carlo Step: %2.0f \n- - - - - - - - - - -\n', RMC)

    for ii = 1:1:length(data1.N)
        fprintf('Matrix Dimension Step: %2.0f \n', data1.N(ii))
        [data1.time_a_method1(RMC,ii)] = a_method1(A{ii},B{ii});
        [data1.time_a_method2(RMC,ii)] = a_method2(A{ii},B{ii});
    end

    fprintf('.\n.\n.\n Elapsed Time: %2.1fs \n', sum([ ... 
        data1.time_a_method1(:);  data1.time_a_method2(:)], 1))
end    

% data1.t_elapsed = sum([data1.time_a_method1(:); data1.time_a_method2(:)], 1);
data1.tEnd = toc(tStart); % pair 2: toc
fprintf('.\n.\n.\nTotal Time: %2.1fs \n', data1.tEnd)


%% (a) Method 1
function [time] = a_method1(A, B)
    tic;
    inv(kron(A,B));
    time = toc;
end


%% (a) Method 2
function [time] = a_method2(A, B)
    tic;
    kron(inv(A),inv(B));
    time = toc;
end