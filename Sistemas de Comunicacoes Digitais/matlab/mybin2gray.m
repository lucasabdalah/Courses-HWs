function gray_const = mybin2gray(bit_seq)
% function gray_const = mybin2gray(M)
    % Computa os simbolos binarios com codificao de Gray, sabendo a quantidade de simbolos (M).
    %
    % SYNTAX: gray_const = mybin2gray(M);
    %
    % INPUTS: 
    %       M : Quantidade de simbolos da constelacao 
    %
    % OUTPUTS:
    %       gray_const : Simbolos binarios com codificao de Gray
    %
    % Referencia: https://www.tutorialspoint.com/conversion-of-binary-to-gray-code
    % 
    %HISTORY:
    % 2021/03/26: - Lucas Abdalah.
    %
    
    % D = de2bi([0:M-1],log2(M),'left-msb') 
    
    K = length(bit_seq);
    g=zeros(1,K);
    b = logical(flip(bit_seq)); %Flip the vector to obtain a LSB 
    n = 0;

    while K > n
        ii=K-n;
        if ii == K
            g(ii) = double(b(ii));
        else
            g(ii) = double(xor(b(ii+1),b(ii)));
        end
        n=n+1;
    end
    
    gray_const = flip(g);

end % function