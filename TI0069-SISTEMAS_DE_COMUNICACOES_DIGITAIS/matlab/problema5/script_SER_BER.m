%% Trabalho de SCD 
% ----------------------------
% Problema 5 - Teorica x AWGN das Modulacoes M-QAM e M-PSK
% script_SER_BER.m
% 2021/03/27 - Lucas Abdalah
%

%% Add necessary folder paths
addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab'
addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab\problema2'
addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab\problema3'
addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab\problema4\parte2\'
addpath 'C:\Users\lukin\Documents\GitHub\Courses-HWs\Sistemas de Comunicacoes Digitais\matlab\problema4\parte3\'

%% General Setup
plot_figures_here = true;
save_figures_here = true;
SER_here = false;
BER_here = true;

if SER_here == true
    fileload = {'Problema2_M_QAM.mat','Problema3_4_QAM.mat','Problema3_16_QAM.mat','Problema3_64_QAM.mat'};

    %% M-QAM
    SER_void_QAM = zeros(6,11);

    load(fileload{1});
    SER_void_QAM(1:3,:) = erro_MQAM;
    close all;
    clearvars -except SER_void_QAM fileload plot_figures_here save_figures_here SER_here

    load(fileload{2});
    SER_void_QAM(4,:) = SER_curve;
    close all;
    clearvars -except SER_void_QAM fileload plot_figures_here save_figures_here SER_here

    load(fileload{3});
    SER_void_QAM(5,:) = SER_curve;
    close all;
    clearvars -except SER_void_QAM fileload plot_figures_here save_figures_here SER_here

    load(fileload{4});
    SER_void_QAM(6,:) = SER_curve;
    close all;
    clearvars -except SER_void_QAM plot_figures_here save_figures_here SER_here

    %% M-PSK
    fileload = {'Problema4_M_PSK.mat','Problema4_4_PSK.mat','Problema4_8_PSK.mat'};

    SER_void_PSK = zeros(4,11);

    load(fileload{1});
    SER_void_PSK(1:2,:) = erro_MPSK;
    close all;
    clearvars -except SER_void_QAM SER_void_PSK fileload plot_figures_here save_figures_here SER_here

    load(fileload{2});
    SER_void_PSK(3,:) = SER_curve;
    close all;
    clearvars -except SER_void_QAM SER_void_PSK fileload plot_figures_here save_figures_here SER_here

    load(fileload{3});
    SER_void_PSK(4,:) = SER_curve;
    close all;
    clearvars -except SER_void_QAM SER_void_PSK plot_figures_here save_figures_here SER_here

    Es_N0 = [0:2:20];
    if plot_figures_here == true
        
        %% M-QAM - SER
        str = {'placeholder','placeholder'};
        h = figure;
        h1 = semilogy_erro(Es_N0,SER_void_QAM(1:3,:),str);
        set(h1,'LineStyle','-'); set(h1,'LineWidth',1);
        hold on
        h2 = semilogy_erro(Es_N0,SER_void_QAM(4:6,:),str);
        set(h2,'Marker','o');
        hold on
        
        %% M-PSK - SER
        str = {'placeholder','placeholder'};
        h1 = semilogy_erro(Es_N0,SER_void_PSK(1:2,:),str);
        set(h1,'LineStyle','-'); set(h1,'LineWidth',1);
        hold on
        h2 = semilogy_erro(Es_N0,SER_void_PSK(3:4,:),str);
        set(h2,'Marker','o');
        
        str = {'$P(e)$ - 4QAM','$P(e)$ - 16QAM','$P(e)$ - 64QAM','SER - 4QAM','SER - 16QAM','SER - 64QAM','$P(e)$ - 4PSK','$P(e)$ - 8PSK','SER - 4PSK','SER - 8PSK'};

        legend(str,'Location','best','interpreter','latex');
        ylabel('$P(e)$ e SER','interpreter','latex');
        title(['$P(e)$ teorica x SER $M$-QMAM e $M$-PSK'],'interpreter','latex');
        ylim([1e-6,1e0]);
        
        pause

        if save_figures_here == true
            %% Save figure 
            saveas(h,['fig/Erro_teoricaxAWGN_Geral'],'pdf');
            
            %% Save two variables, where FILENAME is a variable:
            savefile = 'Problema5_ERRO_Geral.mat';
            save(savefile, 'SER_void_QAM','SER_void_PSK','str');
        end
    end
end

if BER_here == true;
    
    fileload = {'Problema3_4_QAM.mat','Problema3_16_QAM.mat','Problema3_64_QAM.mat'};

    %% M-QAM
    BER_void_QAM = zeros(3,11);

    load(fileload{1});
    BER_void_QAM(1,:) = BER_curve;
    close all;
    clearvars -except BER_void_QAM fileload plot_figures_here save_figures_here BER_here

    load(fileload{2});
    BER_void_QAM(2,:) = BER_curve;
    close all;
    clearvars -except BER_void_QAM fileload plot_figures_here save_figures_here BER_here

    load(fileload{3});
    BER_void_QAM(3,:) = BER_curve;
    close all;
    clearvars -except BER_void_QAM plot_figures_here save_figures_here BER_here

    %% M-PSK
    fileload = {'Problema4_4_PSK.mat','Problema4_8_PSK.mat'};

    BER_void_PSK = zeros(2,11);

    load(fileload{1});
    BER_void_PSK(1,:) = BER_curve;
    close all;
    clearvars -except BER_void_QAM BER_void_PSK fileload plot_figures_here save_figures_here BER_here
    
    load(fileload{2});
    BER_void_PSK(2,:) = BER_curve;
    close all;
    clearvars -except BER_void_QAM BER_void_PSK plot_figures_here save_figures_here BER_here


    Es_N0 = [0:2:20];
    if plot_figures_here == true
        
        %% M-QAM - SER
        str = {'placeholder','placeholder','placeholder'};
        h = figure;
        h1 = semilogy_erro(Es_N0,BER_void_QAM(1:3,:),str);
        set(h1,'Marker','o');
        hold on
        
        %% M-PSK - SER
        str = {'placeholder','placeholder'};
        h1 = semilogy_erro(Es_N0,BER_void_PSK(1:2,:),str);
        set(h1,'Marker','o');
        
        str = {'4QAM','16QAM','64QAM','4PSK','8PSK'};

        legend(str,'Location','best','interpreter','latex');
        ylabel('BER','interpreter','latex');
        title(['BER $M$-QMAM e $M$-PSK'],'interpreter','latex');
        ylim([1e-6,1e0]);
        
        pause

        if save_figures_here == true
            %% Save figure 
            saveas(h,['fig/BER_teoricaxAWGN_Geral'],'pdf');
            
            %% Save two variables, where FILENAME is a variable:
            savefile = 'Problema5_BER_Geral.mat';
            save(savefile, 'BER_void_QAM','BER_void_PSK','str');
        end
    end
end