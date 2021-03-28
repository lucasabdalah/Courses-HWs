function hfig = semilogy_erro(Es_N0,erro_MQAM, str)
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

hfig = semilogy(Es_N0, erro_MQAM,...
'Marker','square',...
'LineWidth', 2.5);
xlabel('$\frac{E_b}{N_0}$ (dB)','interpreter','latex');
ylabel('$P(e)$','interpreter','latex');
legend(str,'Location','best');
xticks(Es_N0);
grid on

end