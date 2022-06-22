%% [TIP8419 - Algebra Linear e Multilinear] - Homework 0 
% by Lucas Abdalah
% ----------------------------
% Problema 1 - Item b: Kronecker Product and Matrix Inversion Cost
% draft.m
% Author: Lucas Abdalah
%


%% General Setup
data2.MC = 100; 
data2.N = 2;
data2.K = 2:2:10;
data2.time_b_method1 = zeros(data2.MC, length(data2.K));
data2.time_b_method2 = zeros(data2.MC, length(data2.K));

for RMC = 1:1:data2.MC
    fprintf('Matrix Dimension: %2.0f --- RMC Step: %2.0f \n', data2.N, RMC)
    for ii = 1:1:length(data2.K)
        fprintf('Number of kron products: %2.0f \n', data2.K(ii))
        [data2.time_b_method1(RMC,ii)] = b_method1(data2.N, data2.K(ii));
        [data2.time_b_method2(RMC,ii)] = b_method2(data2.N, data2.K(ii));
    end
end


%% Save Data
% disp('Saving data...')
% save('hw0_data.mat','-struct','data1');
% save('hw0_data.mat','-struct','data2');
% disp('Saved sucessfully')


%% Gen Figure
close all;
figure()
semilogy(data2.K, mean(data2.time_b_method1, 1), ...
    'Color', 'blue',...        
    'LineStyle', '--',...
    'LineWidth', 1.0,...
    'Marker', 'o',...
    'MarkerFaceColor', 'blue',...
    'MarkerSize', 5);
hold on 
semilogy(data2.K, mean(data2.time_b_method2, 1), ...
    'Color', 'red',...        
    'LineStyle', '-',...
    'LineWidth', 1.0,...
    'Marker', 's',...
    'MarkerFaceColor', 'red',...
    'MarkerSize', 5);
hold off
xticks(data2.K);
xlabel('Matrix Dimension, K')
ylabel('Time (s)')
legend(["Method 1", "Method 2"], 'Location', 'Best') % legend(leg, 'Location', 'Northeastoutside')
legend boxoff
grid on
axis tight


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
