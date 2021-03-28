%%% Simulation of digital modulation techniques
clear 
close all
clc

L=4096;
M= 4;
K= log2(M);
N=L/K;

Es_No=30; % em dB

%%% 4- PAM %%%
Es= 5;
d= sqrt(3*Es/(M^2-1));

% generate bit stream
bits= randi(2,1,L)-1;

% mapping of bits blocks as constellation points
symb= zeros(N,1); k=1;
for i=1:N
    symb(i)=pam4_mapping(bits(k:k+K-1),d);
    k=k+K;
end


%% Conferindo a constela��o transmitida
scatter(symb,zeros(size(symb,1),1),75,[0.25 0.2 0.65],'filled')
xlabel('In-Phase')
ylabel('Quadrature')
title('4-PAM')
% legend('4-PAM constela��o')
grid on
hold on
% 
%% generate the noise term
No= Es*10^(-Es_No/10);
noise= sqrt(No/2)*(randn(N,1)+1i*randn(N,1));

% received signal with noise
rx_signal= symb + noise; 
scatter(real(rx_signal),zeros(size(symb,1),1),15,[0.8 0.6 0.15],'filled')
legend('True Signal','4-PAM constela��o + Noise')
pause;
%% decision (slicer)
% -- Pre-aloca��o de mem�ria
symb_dec=zeros(N,1); % s�mbolo decidido no slicer
bits_dec=zeros(1,L); % bit     decidido no slicer

% - Construa seu slicer
% Aqui  rx_signal  
% [symb_dec,bits_dec] = slicer_4_PAM()

%% C�clulo de BER e SER
% --- Aqui ser�o comparados os s�mbolos (bits) transmitidos originalmente
% "symb" com os decidos symb_dec. Sendo feita a contagem de erro total.

% BER = minha_funcao_BER()
% SER = minha_funcao_SER()

%% Teste do plot P(erro)xEb/N0
%% Calculo da probabilidade do M-QAM

N_size=size(M,2);
Es_N0 = [0:2:20];  % em dB
N_dB = length(Es_N0) ;
N0 = Es'.*(10.^(-Es_N0/10));

% generate the noise term
% noise= sqrt(No/2)*(randn(N,1)+1i*randn(N,1));

erro_MQAM = zeros(N_size,N_dB);

for ii=1:N_size
    for jj=1:N_dB
        erro_MQAM(ii,jj)= Pe_MQAM(M(ii),Es(ii),N0(ii,jj));
    end
end