function output = problem_1
    %problem_1 - Description
    %
    % Syntax: output = myFun(input)
    %
    % Long description
    
    % Item A
    data.MC = 2;
    data.N = 2.^(1:6);
    data.time_a_method1 = zeros(data.MC, length(data.N));
    data.time_a_method2 = zeros(data.MC, length(data.N));

    for ii = 1:1:length(data.N)
        data.A{ii} = complex(randn(data.N(ii), data.N(ii)), randn(data.N(ii),data.N(ii)));
        data.B{ii} = complex(randn(data.N(ii), data.N(ii)), randn(data.N(ii),data.N(ii)));
    end
    
    for RMC = 1:1:data.MC
        
        fprintf('Monte Carlo Step: %2.0e ...\n', RMC)

        for ii = 1:1:length(data.N)
            
            fprintf('Matrix Dimension Step: %2.0e ... \n', ii)

            [data.time_a_method1(RMC,ii), data.C1{RMC,ii}] = a_method1(data.A{ii},data.B{ii});
            [data.time_a_method2(RMC,ii), data.C2{RMC,ii}] = a_method2(data.A{ii},data.B{ii});
        end
    end
    
    
    output = data;

    % Item B

    
end


% function rnd_COMPLEX_MATRIX = gen_COMPLEX_MATRIX(N)
%     rnd_COMPLEX_MATRIX = complex(randn(N, N), randn(N,N));
% end

% Item A
function [time, C] = a_method1(A, B)
    tic;
    C = inv(kron(A,B));
    time = toc;
end


function [time, C] = a_method2(A, B)
    tic;
    C = kron(inv(A),inv(B));
    time = toc;
end

% Item B
function [time, C] = b_method1(K)
    N = 4;
    
    % Generate matrices
    for ii = 1:K
        A{ii} = complex(randn(N, N), randn(N,N));
    end
    
    tic;
    % Produtorio
    for ii = 1:K
        
    end
    time = toc;
end

function [time, C] = b_method2(N)
    tic;

    time = toc;
end