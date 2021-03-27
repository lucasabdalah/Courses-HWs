function constelacao_gray = gray_const(M,verbose)
% function constelacao_gray = gray_const(M)
% Computa a codificao de Gray de um alfabeto binario com log2(M) bits.
%
% SYNTAX: constelacao_gray = gray_const(M);
%
% INPUTS: 
%       M : Quantidade de simbolos a serem representados. (potencia de 2)
%
% OUTPUTS:
%       constelacao_gray : Alfabeto binario em codigo Gray.
%
% EXAMPLE:
%       Exemplo de traducao de codigo binario com 4 bits
%
%       M = 16;
%       D = de2bi([0:M-1],log2(M),'left-msb');
%       G = zeros(M,log2(M));
%       % print a table with both sequences
%       fprintf('Bits\t\tGray\n');
% 
%       for ii=1:M
%           G(ii,:) = gray_const(D(ii,:));
%           disp([num2str(ii), ' ', num2str(D(ii,:)),9,'|',9,num2str(G(ii,:))]);
%       end
%
% Referencia: https://www.tutorialspoint.com/conversion-of-binary-to-gray-code
% 
%HISTORY:
% 2021/03/26: - Lucas Abdalah.
%

D = de2bi([0:M-1],log2(M),'left-msb');
G = zeros(M,log2(M));

for ii=1:M
    G(ii,:) = mybin2gray(D(ii,:));
end
constelacao_gray = G;

% d2 = bi2de(b2,'left-msb')
if verbose == true
    fprintf('\t|Bits|\t\t\t|Gray|\n');
    for ii=1:M
        disp([num2str(ii), 9,'|', num2str(D(ii,:)),'| ->',9,'|', num2str(G(ii,:)),'|']);
    end
end

end % function