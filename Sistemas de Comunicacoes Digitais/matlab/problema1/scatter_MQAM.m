function hfig = scatter_MQAM(r,large_size,received_symbol,write_bits);
% function hfig = scatter_MQAM(M,r);
% Computa a codificao de Gray a partir de uma sequencia de bits de tamanho K.
%
% SYNTAX: hfig = scatter_MQAM(r,large_size,received_symbol);
%
% INPUTS: 
%       r : complex number symbol
%       large_size : to enlarge the circle in the plot
%       received_symbol : to plot a received symbol
% 
% OUTPUTS:
%       hfig : saved figure properties
%
%HISTORY:
% 2021/03/26: - Lucas Abdalah.
%

%% Color definitions for graphics
yellow = [0.9290 0.6940 0.1250];
black  = [0 0 0];
M = size(r,2);
%% Scatter plot
if received_symbol == false
    hfig = scatter(real(r),imag(r), 'filled',...
    'MarkerEdgeColor',yellow,... 
    'MarkerFaceColor',yellow);
    xlabel('In-Phase');
    ylabel('Quadrature');
    title(['Constelacao ',num2str(M),'-QAM']);
    hold on;
    ax = max(real(r))+sqrt(2)/2;
    line([0,0],[-ax,ax], 'Color', black);
    line([-ax,ax],[0,0],'Color', black);
end

if received_symbol == true
    symbol_color=[0.3010 0.7450 0.9330];
    hfig = scatter(real(r),imag(r), 'filled',...
    'Marker','square',...
    'MarkerEdgeColor',symbol_color,... 
    'MarkerFaceColor',symbol_color);
    hfig.SizeData = 20;
end

if large_size == true
    hfig.SizeData = 50;
end

if write_bits == true
% Generating Gray alphabet and constellation's plot.
    for jj = 1:M
        gray_alfabeto = gray_const(M,false);    
        str = {strjoin(string(gray_alfabeto(jj,:)))};
        text(real(r(jj)),imag(r(jj))+0.25* sqrt(2),...
        str,...
        'FontSize', 6,...
        'HorizontalAlignment', 'center');
    end

end


end




%% Desenhar o limite da distancia euclidiana
% Draw a circle radius = d
% center = [real(r) imag(r)];
% radii = sqrt(2)/2;
% viscircles(center,radii,'Color','b');
% axis square
%