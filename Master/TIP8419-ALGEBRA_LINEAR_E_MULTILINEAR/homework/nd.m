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
%       ND.VEC              - Vectorize a matrix.
%       ND.NMSE             - Normalized mean square error (NMSE) of a tensor.
% 
%   MATRIX PRODUCTS
%       ND.HADAMARD_    - Hadamard product with two matrices.
%       ND.KRON_        - Kronnecker product with two matrices.
%       ND.KR_          - Khatri-Rao product with two matrices.
%
%   TENSOR FACTORS ESTIMATION
%       ND.LSKRF        - Least-Squares Khatri-Rao Factorization (LSKRF)
%       ND.LSKRONF      - Least-Squares  Kronecker Product Factorization (LSKRONF)
%       ND.KPSVD        - Kronecker Product Singular Value Decomposition (KPSVD)
%
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

    function y = vec(x)
    % ND.VEC - Vectorize a matrix.
    %   y = vec(x) draws a vector from a given matrix.
    % 
    %   See also.
        y = x(:);
    end

    function [X_nmse, X_nmse_dB] = nmse(X, X_hat)
    % ND.NMSE - Normalized mean square error (NMSE) of a tensor.
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

        if iX == M*N % Verify the input dimensions
            Ahat = complex(zeros(M,jX),0);
            Bhat = complex(zeros(N,jX),0);

            for jj = 1:jX
                [U,S,V] = svd(reshape(X(:,jj), [N M]));
                Ahat(:,jj) = sqrt(S(1,1)).*conj(V(:,1));
                Bhat(:,jj) = sqrt(S(1,1)).*U(:,1);
            end
        else
            error('number of rows of X should be equal to size M*N');
        end
    
    end
    
    
    function [Ahat,Bhat] = lskronf(X, Ma, Na, Mb, Nb)
    % ND.LSKRONF  Least-Squares Kronecker Product Factorization (LSKRONF)
    %   [Ahat,Bhat] = nd.lskronf(X, Ma, Na, Mb, Nb) compute the LSKRONF.
    %
    %   See also.   
        [Mx,Nx] = size(X);
        
        if Ma*Mb == Mx && Na*Nb == Nx % Verify the input dimensions 
            Xhat = complex(zeros(Mb*Nb,Ma*Na),0);    
            X_b = mat2cell(X, repelem(Mx/Ma,Ma), repelem(Nx/Na,Na));
            
            itCol = 1;
            for jj = 1:Na
                for ii = 1:Ma
                    Xhat(:,itCol) = nd.vec(cell2mat(X_b(ii,jj)));
                    itCol = itCol + 1;
                end
            end

            [U,S,V] = svd(Xhat);
            Ahat = reshape(sqrt(S(1,1)).*conj(V(:,1)),[Ma Na]);
            Bhat = reshape(sqrt(S(1,1)).*U(:,1), [Mb Nb]);

        else
            error('size of X(Mx, Nx) should match with Mc=Ma*Mb and Nc=Na*Nb, A(Ma, Na) and B(Mb, Nb)');
        end
    end


    function [U,S,V,rkp] = kpsvd(X, Xstruct)
    % ND.KPSVD  Kronecker Product Singular Value Decomposition (KPSVD)
    %   [U,S,V,rkp] = nd.kpsvd(X, Xstruct) compute the KPSVD.
    %
    %   See also.   
        [Mx,Nx] = size(X);
        
        if Xstruct(1)*Xstruct(3) == Mx && Xstruct(2)*Xstruct(4) == Nx % Verify the input dimensions 
            Xhat = complex(zeros(Xstruct(3)*Xstruct(4),Xstruct(1)*Xstruct(2)),0);    
            X_b = mat2cell(X, repelem(Mx/Xstruct(1),Xstruct(1)), repelem(Nx/Xstruct(2),Xstruct(2)));
            
            itCol = 1;
            for jj = 1:Xstruct(2)
                for ii = 1:Xstruct(1)
                    Xhat(:,itCol) = nd.vec(cell2mat(X_b(ii,jj)));
                    itCol = itCol + 1;
                end
            end
            [U,S,V] = svd(Xhat');
            rkp = rank(S);
        else
            error('size of X(Mx, Nx) should match with Mc=Xstruct(1)*Xstruct(3) and Nc=Xstruct(2)*Xstruct(4), for A(Xstruct(1), Xstruct(2)) and B(Xstruct(3), Xstruct(4))');
        end

    end
    
    

    %% PLACE HOLDER SECTION

end

end