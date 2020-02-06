%Retorna o sinal apos aplicacao da Transformada de Fourrier
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
