%% [TIP8419 - Algebra Linear e Multilinear]
% ----------------------------
% nd.m
% Author: Lucas Abdalah
% 
% ND is a package developped for the Multilinear Algebra Course
% It is a shortcut for N-d array in reference to the homonym library in python
%
% Content
%   MATRIX PRODUCTS
%       ND.HADAMARD_    - Hadamard product with two matrices.
%       ND.KRON_        - Kronnecker product with two matrices.
%       ND.KR_          - Khatri-Rao product with two matrices.
%
%   PLACE HOLDER
%
%   PARAFAC/CP
%

classdef nd

methods(Static)
        
    %% MATRIX PRODUCTS
    function [C, elaspedTime] = hadamard_(A, B)
    % ND.HADAMARD_  Hadamard product with two matrices.
    %   C = nd.hadamard_(A, B) compute the hadamard procuct.
    %
    %   [C, elaspedTime] = nd.hadamard_(A, B) compute the hadamard 
    %   procuct elapsed time.
    %
    %   See also.
        tic;
        
        C = A.*B;
        
        elaspedTime = toc;
    end

    function [C, elaspedTime] = kron_(A, B)
    % ND.KRON_  Kronnecker product with two matrices.
    %   C = nd.kron_(A, B) compute the Kronnecker procuct.
    %
    %   [C, elaspedTime] = nd.kron_(A, B) compute the Kronnecker 
    %   procuct elapsed time.
    %
    %   See also.
        tic;
        
        % [M_rows, N_columns] = size(B);
        % C = repelem(A, M_rows, N_columns).*repmat(B,[size(A)]);
        
        C = kron(A,B);

        elaspedTime = toc;
    end

    function [C, elaspedTime] = kr_(A, B)
    % ND.KR_  Khatri-Rao product with two matrices.
    %   C = nd.kr_(A, B) compute the Khatri-Rao procuct.
    %
    %   [C, elaspedTime] = nd.kr_(A, B) compute the Khatri-Rao 
    %   procuct elapsed time.
    %
    %   See also.
        tic;
        
        N = size(A,2);
        if N == size(B,2)
            P = reshape(A,1,[],N);
            Q = reshape(B,[],1,N);
            C = P.*Q;
            C = reshape(C,[],N);
        else
            error('number of columns should be equal')
        end

        elaspedTime = toc;
    end

    %% PLACE HOLDER SECTION

end

end