%% AP2 de PES
% Questao 3
%
% Q3_script.m
%
% 2021/08/26 - Lucas Abdalah

close all; clearvars; clc; % Clear the matlab ambient

%% To reproduce the same results
rng('default');

%% General setup
L = 20; % Filter Length (Number of coefficients)
delta = 15; % Delay
qam = 16; % 16-QAM
mu = 0.4;
epsilon = 1e-6; 

%% Channel Parameters
H_z=[0.5 1.2 1.5 -1]; % Channel TF 
SNR_dB = 30; % SNR in dB

%% Training vs Decision setup
Ntrain = 500;  % N symbols for 4QAM
Ndecision = 5000; % length of decision-directed with QAM data
N=Ntrain + Ndecision; % total number of symbols
train_data = zeros(1,Ntrain+delta); %training data

s=zeros(1,N);

% Pagina 107 do Proakis (5a ed.)
E_media = 5;
E_simbolo = sqrt(2*E_media); 

%% SNR's for both cases
SNR_train = sqrt(norm(H_z)^2/10^(SNR_dB/10)); % during training
SNR_decision = sqrt((E_simbolo)^2*norm(H_z)^2/10^(SNR_dB/10)); % during decision-directed

%% 
v=zeros(1,N);
v(1:Ntrain)=(SNR_train^2/2)*(randn(1,Ntrain)+1i*randn(1,Ntrain));
v(Ntrain:N)=(SNR_decision^2/2)*(randn(1,N-Ntrain+1)+1i*randn(1,N-Ntrain+1));

% 4QAM sequence for trainining
s(1:Ntrain)=(sign(randn(1,Ntrain))+1i*sign(randn(1,Ntrain)))/(sqrt(2)); % 
train_data(delta+1:Ntrain+delta)=s(1:Ntrain); 
   
for i=1:Ndecision
 xint=randi([0,(qam-1)],1);
 xcoor=qammod(xint,qam);% Matlab Mod function for decision 
 s(Ntrain+i)=xcoor; 
end 
 
% Pass the data through the channel
y=filter(H_z,1,s); 
r=y+v; 

%% LMS algorithm  
w  = zeros(L,1); %
u  = zeros(1,L); 
e = zeros(1,N); % error vector
count_errors=0;

% Training mode
for i = 1:Ntrain+delta
   u  = [r(i) u(1:L-1)];
   hat_s(i) = u*w;
   d(i) = train_data(i); % training
   e(i) = d(i)- hat_s(i);
   w = w + mu*e(i)*u'/(norm(u)^2+epsilon);
end

% Decision-directed mode
for i = Ntrain+delta+1:N
   u  = [r(i) u(1:L-1)];
   hat_s(i) = u*w;
   
   check_s(i)=slicer_16QAM(hat_s(i));
   
   d(i) = check_s(i); % decision_directed
   e(i) = d(i)- hat_s(i);
   w = w + mu*e(i)*u'/(norm(u)^2+epsilon);
   
   if (check_s(i) ~= s(i-delta))
     count_errors = count_errors + 1;
   end
   
end  

figure
subplot(221); 
plot(real(train_data(delta+1:Ntrain)),imag(train_data(delta+1:Ntrain)),'.');
grid on; 
str = ['Training']; title(str,'interpreter','latex');

subplot(222); 
plot(real(s(Ntrain+delta+1:N)),imag(s(Ntrain+delta+1:N)),'.');
grid on; 
str = ['Transmitted']; title(str,'interpreter','latex');

subplot(223); 
plot(real(r(Ntrain+delta+1:N)),imag(r(Ntrain+delta+1:N)),'.');grid; 
grid on;  
str = ['Received']; title(str,'interpreter','latex');

subplot(224); 
plot(real(hat_s(Ntrain+delta+1:N)),imag(hat_s(Ntrain+delta+1:N)),'.');grid; 
grid on;
str = ['Output']; title(str,'interpreter','latex');
