%% [TIP8419 - Algebra Linear e Multilinear]
% ----------------------------
% Author: Lucas Abdalah
% nd.m
% 
% ND is a package developped for the Multilinear Algebra Course
% It is a shortcut for N-d array in reference to the homonym library in python
%
% Content
%   MATRIX FUNCTIONS
%       ND.RANDN_COMPLEX    - Complex-valued matrix from normal distribution.
% 
%   MATRIX PRODUCTS
%       ND.HADAMARD_    - Hadamard product with two matrices.
%       ND.KRON_        - Kronnecker product with two matrices.
%       ND.KR_          - Khatri-Rao product with two matrices.
%
%   TENSOR FACTORS ESTIMATION
%       ND.LSKRF        - Least-Squares Khatri-Rao Factorization (LSKRF)
%
%   PARAFAC/CP
%

classdef nd

methods(Static)
    
    %% MATRIX FUNCTIONS
    function C = randn_complex(M, N)
    % ND.RANDN_COMPLEX - Complex-valued matrix from normal distribution.
    %   C = nd.randn_complex(M,N) draws a complex-valued matrix from normal 
    %       distribution.
    % 
    %   See also.
        C = complex(randn(M,N), randn(M,N));
    end

    function [X_nmse, X_nmse_dB] = nmse(X, X_hat)
    % ND.nmse - Normalized mean square error (NMSE) of a tensor.
    %   [X_nmse, X_nmse_dB] = nd.nmse(X, X_hat) compute the NMSE of two arrays
    % 
    %   See also.
        X_nmse = frob(X - X_hat)^2/(frob(X)^2);
        X_nmse_dB = 20*log10(X_nmse);
    end


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


    %% TENSOR FACTORS ESTIMATION
    function [Ahat,Bhat] = lskrf(X, M, N)
    % ND.LSKRF  Least-Squares Khatri-Rao Factorization (LSKRF)
    %   [Ahat,Bhat] = nd.lskrf(X, M, N) compute the LSKRF.
    %
    %   See also.
        [iX, jX] = size(X);

        if iX == M*N
            Ahat = complex(zeros(M,jX),0);
            Bhat = complex(zeros(N,jX),0);

            for jj = 1:jX
                Xp = reshape(X(:,jj), [N M]);
                [U,S,V] = svd(Xp);
                Ahat(:,jj) = sqrt(S(1,1)).*conj(V(:,1));
                Bhat(:,jj) = sqrt(S(1,1)).*U(:,1);
            end
        else
            error('number of rows of X should be equal to size M*N');
        end
    
    end
    

    %% PLACE HOLDER SECTION

end

end

% targetSignals = []; %Default
% for ii = 1:2:numel(varargin)
%     switch lower(varargin{ii})
%         case 'assigntovariables'
%             assignToVariables = varargin{ii+1};
%         case 'targetsignals'
%             targetSignals = varargin{ii+1};
%         otherwise
%             error('EDFREAD: Unrecognized parameter-value pair specified. Valid values are ''assignToVariables'' and ''targetSignals''.')
%     end
% end
