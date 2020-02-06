%PROBLEMAS 1 e 2

% (PARTE 1)  Separacao do Male1_speech 
% Aqui carregamos os vetores
load('Mic1.mat'); 
load('Mic2.mat');
load('Male1_speech.mat');

% Fazendo a transposicao e "corte" do vetor y 
y = y';
y = y(1:length(Mic1));

%A partir de agora utilizamos os calculos demonstrados anteriormente
%Mic1
for i = 1:1:length(Mic1)    
  c(i) = y(i)*Mic1(i);
end;
cory1 = mean(c);
k1 = cory1/var(y);
M1 = Mic1 - k1*y;
%Aqui temos o Mic1 sem a voz do Locutor

%Mic2
for i = 1:1:length(Mic2)
  d(i) = y(i)*Mic2(i);
end;
cory2 = mean(d);
k2 = cory2/var(y);
M2 = Mic2 - k2*y;
%Aqui temos o Mic2 sem a voz do Locutor

% (PARTE 2)  Aplicacoes Algebricas
%Calculo dos coeficientes de correlacao
e11 = corr(M1,M1);
e12 = corr(M1,M2);
e21 = corr(M2,M1);
e22 = corr(M2,M2);
%Aqui construimos a matriz dos coeficientes de Correlacao

%entre os microfones (Sem o sinal de Y)
M = [e11 e12; e21 e22];

%Fazemos a opecao da raiz da Matriz M
M = M^(-0.5);

%Agora utilizamos os coeficientes das matrizes obtidas
X1 = (M(1,1))*M1 + (M(1,2))*M2;
X2 = (M(2,1))*M1 + (M(2,2))*M2; 


%PROBLEMA 3
acy = xcorr(y);
acM1 = xcorr(X1);
acM2 = xcorr(M2);

%Plot das autocorrelacoes
%y
plot(acy);
title('Autorrelacao de y');
xlabel('Amostras');
%X1
plot(X1);
title('Autorrelacao de X1');
xlabel('Amostras');
%X2
plot(X2);
title('Autorrelacao de X2');
xlabel('Amostras');


%PROBLEMA 4
%DEP e Utilizacao de Fourrier
%Primeiramente descrevemos a criacao da funcao da Transformada
%Sendo passsado o vetor da xcorr e a frequencia de amostragem
function [S, frequency] = transformadaDEP(s, fs)
normal = length(s);
aux = 0:normal-1;
T = normal/fs;
frequency = aux/T;
S = fftn(s)/normal;
fc = ceil(normal/2);
S = S(1:fc);
%Feita a plotagem da DEP e a descricao dos eixos
figure();
plot(frequency(1:fc),abs(S));
title('DEP');
xlabel('Frequência em Hertz');
ylabel('Amplitude');

%Depois utilizamos a funcao que ja calcula e plota a DEP
 transformadaDEP(y, Fs);
 transformadaDEP(X1, Fs)
 transformadaDEP(X2, Fs)