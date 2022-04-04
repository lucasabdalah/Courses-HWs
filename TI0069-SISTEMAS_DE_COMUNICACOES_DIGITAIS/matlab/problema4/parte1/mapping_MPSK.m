% Problema 1 - Parte 4
function symb = mapping_MPSK(bits)
% function  symb= mapping_MPSK(bits,M,E_g)
% Computa o mapa da constelacao M-PSK (Gray) a partir de uma sequencia de bits de tamanho M e energia de pulso/modulacao E_g.
%
% SYNTAX: symb= mapping_MPSK(bits,M,E_g);
%
% INPUTS: 
%       bits : Sequencia de bits, Simbolo binarios
%       M : Quantidade de simbolos da constelacao 
%       E_g : Energia de pulso/modulacao E_g.
% OUTPUTS:
%       symb : Simbolo da constelacao
%
% EXAMPLE:
%       bits = [1,0,1];
%       symb = mapping_MPSK(bits);
%
%
%HISTORY:
% 2021/03/26: - Lucas Abdalah.
%

%% A priori info
% --------------------------------------
K = length(bits); M = 2^K; 
E_g = 1; energia_media_PSK = energia_MPSK(M,E_g); 
d = d_MPSK(M,energia_media_PSK);
% --------------------------------------
%% Creating a Gray alphabet for matching the M-PSK constellation
const_PSK = const_MPSK(M,energia_media_PSK);            % Constellation M-PSK
gray_alfabeto = gray_const(M,false);    % Gray alphabet

% To compare gray_alfabeto with the entry bits, we can use xor(,) and get the line where all logical elements are equal.
% get the xor_out' and the position equal to ->size(bits) stands for the correct element
get_position = find(sum(xor(bits,gray_alfabeto)')==0); % Condition to find the equivalent on the Gray dicionary
symb = const_PSK(get_position);

end % function