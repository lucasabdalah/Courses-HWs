%% Item (b)
% N = 200; J = 30;
% A = randn(N,N);
% B = randn(J,J);
% kron_dim(A,B);


%% General Setup
data2.MC = 1; 
data2.N = 4;
data2.K = [2 4 6 8 10];
data2.time_b_method1 = zeros(data2.MC, length(data2.K));
data2.time_b_method2 = zeros(data2.MC, length(data2.K));


[data2.time_b_method1] = b_method1(data2.N,8);


%% Save Data
% disp('Saving data...')
% save('hw0_data.mat','-struct','data1');
% save('hw0_data.mat','-struct','data2');
% disp('Saved sucessfully')


%% (b) Method 1
function [time] = b_method1(N, K)
    fprintf("Start (b) method 1\n")
    A = cell(K,1);
    
    % fprintf("Matrix Dimensions: %dX%d \nN of elements: %d \nMemory use: %2.0f Mb \n", ...
    %     dim(1), dim(2), elements, memory/1024);

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
            kron_dim(C, A{ii})
            C = kron(C, A{ii});
        end
    end
    C = inv(C);
    time = toc;
end


%% (b) Method 2
