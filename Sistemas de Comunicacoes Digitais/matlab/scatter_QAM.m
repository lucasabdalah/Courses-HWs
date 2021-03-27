function hfig = scatter_QAM(const_QAM,large_size);
% function hfig = scatter_QAM(M,const_QAM);
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

%% Color definitions for graphics
yellow = [0.9290 0.6940 0.1250];
black  = [0 0 0];
M = size(const_QAM,2);
%% Scatter plot
hfig = scatter(real(const_QAM),imag(const_QAM), 'filled',...
'MarkerEdgeColor',yellow,... 
'MarkerFaceColor',yellow);
xlabel('In-Phase');
ylabel('Quadrature');
title(['Constelacao ',num2str(M),'-QAM']);
hold on;
ax = max(real(const_QAM))+sqrt(2)/2;
line([0,0],[-ax,ax], 'Color', black);
line([-ax,ax],[0,0],'Color', black);

if large_size == true
    hfig.SizeData = 200;
end

% end