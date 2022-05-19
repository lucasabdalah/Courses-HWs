%% Homework 0 - Multilinear Algebra
% ----------------------------
% Problema 1 - Item a: Kronecker Product and Matrix Inversion Cost
% hw0_problem1.m
% Author: Lucas Abdalah
%

%% General Setup
data.MC = 100;
data.N = 2.^(1:6);
data.time_a_method1 = zeros(data.MC, length(data.N));
data.time_a_method2 = zeros(data.MC, length(data.N));
A = cell(data.MC, length(data.N));
B = cell(data.MC, length(data.N));

tStart = tic;

for ii = 1:1:length(data.N)
 
    A{ii} = complex(randn(data.N(ii), data.N(ii)), randn(data.N(ii),data.N(ii)));
    B{ii} = complex(randn(data.N(ii), data.N(ii)), randn(data.N(ii),data.N(ii)));

end


for RMC = 1:1:data.MC
    
    fprintf('.\n.\n.\nMonte Carlo Step: %2.0f \n- - - - - - - - - - -\n', RMC)

    for ii = 1:1:length(data.N)
        
        fprintf('Matrix Dimension Step: %2.0f \n', data.N(ii))

        [data.time_a_method1(RMC,ii)] = a_method1(A{ii},B{ii});
        [data.time_a_method2(RMC,ii)] = a_method2(A{ii},B{ii});
    end

    fprintf('.\n.\n.\n Elapsed Time: %2.1fs \n', sum([ ... 
        data.time_a_method1(:);  data.time_a_method2(:)], 1))
    
end    

% data.t_elapsed = sum([data.time_a_method1(:); data.time_a_method2(:)], 1);
data.tEnd = toc(tStart); % pair 2: toc
fprintf('.\n.\n.\nTotal Time: %2.1fs \n', data.tEnd)


%% Save Data
% disp('Saving data...')
% save('hw0_data.mat','-struct','data');
% disp('Saved sucessfully')

% Method 1
function [time] = a_method1(A, B)
    tic;
    inv(kron(A,B));
    time = toc;
end

% Method 2
function [time] = a_method2(A, B)
    tic;
    kron(inv(A),inv(B));
    time = toc;
end