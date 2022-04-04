%% Trabalho de SCD 
% ----------------------------
% Problema 3 - Teorica x AWGN da Modulacao M-QAM
% script_teoricaxAWGN.m
% 2021/03/27 - Lucas Abdalah
%

%% Add necessary folder paths
addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab'
addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab\problema2'

%% General Setup
plot_figures_here = true;
save_figures_here = true;

fileload = {'Problema2_M_QAM.mat','Problema3_4_QAM.mat','Problema3_16_QAM.mat','Problema3_64_QAM.mat'};

SER_void = zeros(6,11);

load(fileload{1});
SER_void(1:3,:) = erro_MQAM;
close all;
clearvars -except SER_void fileload plot_figures_here save_figures_here

load(fileload{2});
SER_void(4,:) = SER_curve;
close all;
clearvars -except SER_void fileload plot_figures_here save_figures_here

load(fileload{3});
SER_void(5,:) = SER_curve;
close all;
clearvars -except SER_void fileload plot_figures_here save_figures_here

load(fileload{4});
SER_void(6,:) = SER_curve;
close all;
clearvars -except SER_void plot_figures_here save_figures_here

Es_N0 = [0:2:20];
if plot_figures_here == true
    str = {'placeholder','placeholder','placeholder'};
    h = figure;
    h1 = semilogy_erro(Es_N0,SER_void(1:3,:),str);
    set(h1,'LineStyle','-'); set(h1,'LineWidth',1);
    hold on
    h2 = semilogy_erro(Es_N0,SER_void(3:6,:),str);
    set(h2,'Marker','o');
    
    str = {'$P(e)$ - 4QAM','$P(e)$ - 16QAM','$P(e)$ - 64QAM','SER - 4QAM','SER - 16QAM','SER - 64QAM'};
    legend(str,'Location','best','interpreter','latex');
    
    ylabel('$P(e)$ e SER','interpreter','latex');
    title(['$P(e)$ teorica x SER $M$-QAM'],'interpreter','latex');
    ylim([1e-6,1e0]);
    
    pause

    if save_figures_here == true
        %% Save figure 
        saveas(h,['fig/Erro_teoricaxAWGN_MQAM'],'pdf');
        
        %% Save two variables, where FILENAME is a variable:
        savefile = 'Problema3_ERRO_M_QAM.mat';
        save(savefile, 'SER_void', 'str');
    end
end
