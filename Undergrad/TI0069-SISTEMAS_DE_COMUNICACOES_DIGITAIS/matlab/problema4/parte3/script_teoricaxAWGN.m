%% Trabalho de SCD 
% ----------------------------
% Problema 4 - Teorica x AWGN da Modulacao M-PSK
% script_teoricaxAWGN.m
% 2021/03/27 - Lucas Abdalah
%

%% Add necessary folder paths
addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab'
addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab\problema4\parte2\'

%% General Setup
plot_figures_here = true;
save_figures_here = true;

fileload = {'Problema4_M_PSK.mat','Problema4_4_PSK.mat','Problema4_8_PSK.mat'};

SER_void = zeros(4,11);

load(fileload{1});
SER_void(1:2,:) = erro_MPSK;
close all;
clearvars -except SER_void fileload plot_figures_here save_figures_here

load(fileload{2});
SER_void(3,:) = SER_curve;
close all;
clearvars -except SER_void fileload plot_figures_here save_figures_here

load(fileload{3});
SER_void(4,:) = SER_curve;
close all;
clearvars -except SER_void plot_figures_here save_figures_here

Es_N0 = [0:2:20];
if plot_figures_here == true
    str = {'placeholder','placeholder'};
    h = figure;
    h1 = semilogy_erro(Es_N0,SER_void(1:2,:),str);
    set(h1,'LineStyle','-'); set(h1,'LineWidth',1);
    hold on
    h2 = semilogy_erro(Es_N0,SER_void(3:4,:),str);
    set(h2,'Marker','o');
    
    str = {'$P(e)$ - 4PSK','$P(e)$ - 8PSK','SER - 4PSK','SER - 8PSK'};
    legend(str,'Location','best','interpreter','latex');
    
    ylabel('$P(e)$ e SER','interpreter','latex');
    title(['$P(e)$ teorica x SER $M$-PSK'],'interpreter','latex');
    ylim([1e-6,1e0]);
    
    pause

    if save_figures_here == true
        %% Save figure 
        saveas(h,['fig/Erro_teoricaxAWGN_MPSK'],'pdf');
        
        %% Save two variables, where FILENAME is a variable:
        savefile = 'Problema3_ERRO_M_PSK.mat';
        save(savefile, 'SER_void', 'str');
    end
end
