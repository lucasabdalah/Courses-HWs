function hfig = semilogy_erro(Es_N0,erro, str)
    % function hfig = semilogy_erro(Es_N0,erro, str);
    % Computa a codificao de Gray a partir de uma sequencia de bits de tamanho K.
    %
    % SYNTAX: hfig = semilogy_erro(Es_N0,erro, str)
    % 
    % OUTPUTS:
    %       hfig : saved figure properties
    %
    %HISTORY:
    % 2021/03/26: - Lucas Abdalah.
    %

hfig = semilogy(Es_N0, erro,...
                'Marker','square',...
                'LineWidth', 1.5,...
                'LineStyle', ':');
xlabel('$\frac{E_b}{N_0}$ (dB)','interpreter','latex');
ylabel('$P(e)$','interpreter','latex');
legend(str,'Location','best');
xticks(Es_N0);
grid on

end