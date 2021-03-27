%%% Simulation of digital modulation techniques
clearvars;clc

L=8;
M= 4;
K= log2(M);
N=L/K;

Es_No=0; % em dB

%%% 4- PAM %%%
Es= 1;
d= sqrt(3*Es/(M^2-1));

pause

% generate bit stream
bits= randi(2,1,L)-1;

% mapping of bits blocks as constellation points
symb= zeros(N,1); k=1;
for i=1:N
    symb(i)=pam4_mapping(bits(k:k+K-1),d);
    k=k+K;
end

pause 

% generate the noise term
No= Es*10^(-Es_No/10);
noise= sqrt(No/2)*(randn(N,1)+1i*randn(N,1));

% received signal with noise
rx_signal= symb + noise; 

% decision (slicer)
symb_dec=zeros(N,1);
bits_dec=zeros(1,L);
k=1;
for i=1:N
    [symb_dec(i),bits_dec(k:k+K-1)]= pam4_demapping(real(rx_signal(i)),d);
    k=k+K;
end

diff_vec=symb-symb_dec;
SER= sum(diff_vec~=0)/N

BER= sum(xor(bits,bits_dec))/L
