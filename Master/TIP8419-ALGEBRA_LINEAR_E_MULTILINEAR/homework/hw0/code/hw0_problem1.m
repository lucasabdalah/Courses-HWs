%% Homework 0 - Multilinear Algebra
% ----------------------------
% Problema 1 - Item a: Kronecker Product and Matrix Inversion Cost
% hw0_problem1.m
% Author: Lucas Abdalah
%

workspace = 'b';

if workspace == 'a'

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

end 

if workspace == 'b'
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

end

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


%% kron dimension
function [dim, elements, memory] = kron_dim(A, B)
    % Check if it's real
    if isreal(A) && isreal(B) == true
        memory_elem = 8; % bytes
    else
        memory_elem = 16; % bytes
    end
    
    dim = size(A).*size(B); % get_limit = memory % get_limit.MaxPossibleArrayBytes % 13999538176 bytes
    elements = prod(dim(:));
    memory = (memory_elem*elements);
    
    if memory < 1025
        memory_use = sprintf('%2.0f b\n', memory);
    elseif memory < 1024^2+1
        memory_use = sprintf('%2.0f Kb\n', memory/1024);
    elseif memory < 1024^3+1
        memory_use = sprintf('%2.0f Mb\n', memory/1024^2);
    else
        memory_use = sprintf('%2.0f Gb\n', memory/1024^3);
    end
 
    fprintf("Matrix Dimensions: %dX%d \nN of elements: %d \nMemory use: %s ", ...
        dim(1), dim(2), elements, memory_use);
end