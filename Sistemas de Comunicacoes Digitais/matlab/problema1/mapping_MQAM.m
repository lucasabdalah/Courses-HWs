function symb = mapping_MQAM(bits)
% function  symb= mapping_MQAM(bits,M,E_g)
% Computa o mapa da constelacao M-QAM (Gray) a partir de uma sequencia de bits de tamanho M e energia de pulso/modulacao E_g.
%
% SYNTAX: symb= mapping_MQAM(bits,M,E_g);
%
% INPUTS: 
%       bits : Sequencia de bits, Simbolo binarios
%       M : Quantidade de simbolos da constelacao 
%       E_g : Energia de pulso/modulacao E_g.
% OUTPUTS:
%       symb : Simbolos binarios com codificao de Gray
% EXAMPLE:
%       bits = [1,0,1,1];
%       symb = mapping_MQAM(bits);
%
%
%HISTORY:
% 2021/03/26: - Lucas Abdalah.
%

%% A priori info
% --------------------------------------
K = length(bits); M = 2^K; 
E_g = 1; energia_media_QAM = energia_MQAM(M,E_g); 
d = d_MQAM(M,energia_media_QAM);
% --------------------------------------
%% Creating a Gray alphabet for matching the M-QAM constellation
const_QAM = const_MQAM(M,d);        % Constellation M-QAM
gray_alfabeto = gray_const(M,true);    % Gray alphabet

% To compare gray_alfabeto with the entry bits, we can use xor(,) and get the line where all logical elements are equal.
% get the xor_out' and the position equal to ->size(bits) stands for the correct element
get_position = find(sum(xor(bits,gray_alfabeto)')==0); % Condition to find the equivalent on the Gray dicionary
symb = const_QAM(get_position);

end % function