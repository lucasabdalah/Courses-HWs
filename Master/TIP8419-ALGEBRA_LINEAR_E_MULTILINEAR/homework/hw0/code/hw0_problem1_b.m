%% [TIP8419 - Algebra Linear e Multilinear] - Homework 0 
% by Lucas Abdalah
% ----------------------------
% Problema 1 - Item (b): Kronecker Product and Matrix Inversion Cost
% hw0_problem1_b.m
% Author: Lucas Abdalah
%


%% General Setup
data2.MC = 200; 
data2.N = 2;
data2.K = 2:2:10;
data2.time_b_method1 = zeros(data2.MC, length(data2.K));
data2.time_b_method2 = zeros(data2.MC, length(data2.K));

tStart = tic;


%% Matrix Computation
fprintf('Matrix Dimension: %2.0f \n.\n', data2.N)
for RMC = 1:1:data2.MC
    fprintf('Monte Carlo Step: %2.0f \n- - - - - - - - - - -\n', RMC)

    for ii = 1:1:length(data2.K)
        % fprintf('Number of kron products: %2.0f \n', data2.K(ii))
        [data2.time_b_method1(RMC,ii)] = b_method1(data2.N, data2.K(ii));
        [data2.time_b_method2(RMC,ii)] = b_method2(data2.N, data2.K(ii));
    end
    
    fprintf('\n Elapsed Time: %2.1fs \n', sum([ ... 
        data2.time_b_method1(:); data2.time_b_method2(:)], 1))
end

fprintf('.\n.\n.\n Method 1 (Elapsed Time): %2.1fs \n', sum(...
    data2.time_b_method1(:), 1) )

fprintf('.\n.\n.\n Method 2 (Elapsed Time): %2.1fs \n', sum(... 
    data2.time_b_method2(:), 1))

data2.tEnd = toc(tStart); % pair 2: toc
fprintf('.\n.\n.\nTotal Time: %2.1fs \n.\n.\n', data2.tEnd)


%% Save Data
disp('Saving data...')
save('hw0_data2.mat','-struct','data2');
disp('Saved sucessfully')


%% (b) Method 1
function [time] = b_method1(N, K)
    % fprintf("Start (b) method 1\n")
    A = cell(K,1);
    
    % Generate matrices
    for ii = 1:1:K 
        A{ii} = complex(randn(N, N), randn(N,N));
    end 

    tic;
    % Product
    for ii = 1:1:K
        if ii == 1
            C = A{ii};
        else
            % ram_use.kron_dim(C, A{ii})
            C = kron(C, A{ii});
        end
    end
    C = inv(C);
    time = toc;
end


%% (b) Method 2
function [time] = b_method2(N, K)
    % fprintf("Start (b) method 2\n")
    A = cell(K,1);
    
    % Generate matrices
    for ii = 1:1:K 
        A{ii} = complex(randn(N, N), randn(N,N));
    end 

    tic;
    % Product
    for ii = 1:1:K
        if ii == 1
            C = inv(A{ii});
        else
            % ram_use.kron_dim(C, A{ii})
            C = kron(C, inv(A{ii}));
        end
    end
    time = toc;
end
