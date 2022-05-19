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
        
        fprintf('Monte Carlo Step: %2f ...\n', RMC)

        for ii = 1:1:length(data.N)
            
            fprintf('Matrix Dimension Step: %2f ... \n', ii)

            [data.time_a_method1(RMC,ii)] = a_method1(data.A{ii},data.B{ii});
            [data.time_a_method2(RMC,ii)] = a_method2(data.A{ii},data.B{ii});
        end
    end    
    
    output = data;

    % Item B

    % data2.MC = 2;
    % data2.N = 2.^(1:6);
    % data2.time_a_method1 = zeros(data.MC, length(data.N));
    % data2.time_a_method2 = zeros(data.MC, length(data.N));



%     K = 2:2:10;
%     N = length(K);
%     time_cost = zeros(MC, N);

%     for ii = 1:1:N
%         for jj = 1:1:K(ii)
%             time_cost(ii) = b_method1(K);
%             disp(jj)
%         end
%     end

%     output = time_cost;
    
% end

% Item A
% function [time] = a_method1(A, B)
%     tic;
%     inv(kron(A,B));
%     time = toc;
% end

% function [time] = a_method2(A, B)
%     tic;
%     kron(inv(A),inv(B));
%     time = toc;
% end

% % Item B
function [time] = b_method1(K)
    N = 4;
    C = ones(4,4);
    A = cell(K,1);
    % Generate matrices
    for ii = 1:K
        A{ii} = complex(randn(N, N), randn(N,N));
    end
    
    tic;
    % Produtorio
    for ii = 1:1:K
        C = C*A{ii};
    end
    time = toc;
end

% function [time, C] = b_method2(N)
%     tic;

%     time = toc;
% end

% function complex_matrix = randn_complex_matrix(N)
%     complex_matrix = complex(randn(N, N), randn(N,N));
% end