function gray_seq = mybin2gray(bit_seq)
% function gray_seq = mybin2gray(M)
% Computa a codificao de Gray a partir de uma sequencia de bits de tamanho K.
%
% SYNTAX: gray_seq = mybin2gray(M);
%
% INPUTS: 
%       M : Quantidade de simbolos da constelacao 
%
% OUTPUTS:
%       gray_seq : Simbolos binarios com codificao de Gray
%
% EXAMPLE:
%       Exemplo de traducao de uma sequencia binaria de 4 bits
%
%       bit_len = 4;
%       b = randi([0 1],bit_len,1)';
%       g = mybin2gray(b);
%       disp(['b ==> ', num2str(b)])
%       disp(['g ==> ', num2str(g)])
%
% Referencia: https://www.tutorialspoint.com/conversion-of-binary-to-gray-code
% 
%HISTORY:
% 2021/03/26: - Lucas Abdalah.
%
        
K = length(bit_seq); % Tamanho da Sequencia de Bits
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

gray_seq = flip(g);

end % function