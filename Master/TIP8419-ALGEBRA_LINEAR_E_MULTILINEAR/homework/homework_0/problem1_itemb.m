%% Homework 0 - Multilinear Algebra
% ----------------------------
% Problema 1 - Cost: Kronecker Product and Matrix Inversion 
% problem1_itemb.m
% Author: Lucas Abdalah
%

% Problem 1 - For randomly generated $\mathbf{A} \in \mathbb{C}^{N \times N}$ and $\mathbf{B} \in \mathbb{C}^{N \times N}$ 


N = 4;
K = [2 4 6 8 10];
time1 = zeros(length(K),1);
time2 = zeros(length(K),1);

[time, ~] = b_method1(7)
% K = 2:2:10;
% N = length(K);
% time_cost = zeros(MC, N);

% for ii = 1:1:N
%     for jj = 1:1:K(ii)
%         time_cost(ii) = b_method1(K);
%         disp(jj)
%     end
% end

function [time, C] = b_method1(K)
    N = 4;
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
            C = kron(C, A{ii});
        end
    end
    C = inv(C);
    time = toc;
    disp(size(C))
end


% x = 1:10;
% n = length(x);
% avg = mymean(x,n)
% med = mymedian(x,n)

% function a = mymean(v,n)
% % MYMEAN Local function that calculates mean of array.

%     a = sum(v)/n;
% end

% function m = mymedian(v,n)
% % MYMEDIAN Local function that calculates median of array.

%     w = sort(v);
%     if rem(n,2) == 1
%         m = w((n + 1)/2);
%     else
%         m = (w(n/2) + w(n/2 + 1))/2;
%     end
% end