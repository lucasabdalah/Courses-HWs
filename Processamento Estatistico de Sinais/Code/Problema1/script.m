%% AP2 de PES
% 2021/08/23 - Lucas Abdalah
%
%close all; clearvars; clc; % Clear the matlab ambient

%%
%
% https://www.mathworks.com/matlabcentral/fileexchange/3649-lms-algorithm-demo
% https://www.mathworks.com/matlabcentral/fileexchange/35670-lms-algorithm-implementation
% https://www.mathworks.com/matlabcentral/fileexchange/70139-lms-algorithm
%
%
%{
https://www.mathworks.com/help/dsp/ref/dsp.lmsfilter-system-object.html#d123e264963
https://www.mathworks.com/help/dsp/ug/system-identification-fir-filter-using-lms-algorithm.html
https://www.mathworks.com/help/dsp/ug/enhance-a-signal-using-lms-and-normalized-lms-algorithms.html
%}
%% 
% disp('testtesttesttesttesttesttesttesttest')

clc;
close all;
clear all;
t=0.001:0.001:1;
d=2*sin(2*pi*50*t);
figure;M=25;
plot(d);N=numel(d);
x=d(1:N)+.9*randn(1,N);
w=zeros(1,M);wi=zeros(1,M);
e=[];mu=1e-4;% u is taken small to ensure that the algorithm converges
for i=M:N
e(i)=d(i)-wi*x(i:-1:i-M+1)';
wi = wi+2*mu*e(i)*x(i:-1:i-M+1);
end
y=zeros(N,1);
for i=M:N
  j=x(i:-1:i-M+1);
  y(i)=((wi)*(j)');
end
ee=y'-d;
subplot(4,1,1),plot(d);ylim([-5 5]);title('desired signal');
subplot(4,1,2),plot(x);title('signal corrupted with noise');
subplot(4,1,3),plot(y);title('estimated signal');
subplot(4,1,4),plot(ee);ylim([-5 5]);title('error signal');
figure(); 
error_1_dB = 10*log10(ee.^2);
subplot(2,1,1);plot(ee.^2);title('error signal ee^2');
subplot(2,1,2);plot(error_1_dB);title('MSE dB');
