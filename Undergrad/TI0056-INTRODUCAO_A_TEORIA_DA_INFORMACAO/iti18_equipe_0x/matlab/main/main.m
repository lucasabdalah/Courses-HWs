lena = imread('lenaTest2.jpg'); %adquire imagem Lena como matriz
tam = size(lena); %obter tamanho da matriz da imagem
lena_v = double(reshape(lena,tam(1)*tam(2),1)); %reordena a matriz como um vetor
lena_b = de2bi(lena_v,8); %converte de decimal para binário
tam_b = size(lena_b); %obtem tamanho da nova matriz com dados em binário
pe = 0:.25:1; %vetor com probabilidades de erro de canal pe
bsc_im = zeros([tam_b length(pe)]); %pré alocando array com imagens pós passagem pelo canal
for p = 1:length(pe) %laço para passagens do canal aumentando gradativamente o a probabilidade de eerro de transmissão
    bsc_im(:,:,p) = bsc(lena_b,pe(p)); %função que simula canal bsc com para um erro pe
end
for p = 1:length(pe) %laço para vizualir resultados
    im = uint8(reshape(bi2de(bsc_im(:,:,p)),tam)); %transformando imagens para seu formato otiginal
    figure %plotando resultados
    imshow(im) 
    xlabel(['$p_e = ' num2str(pe(p)*100)  '\%$' ],'Interpreter','Latex','FontSize',18) %legenda
    saveas(gcf,[num2str(pe(p)*100) 'errolena.pdf']) %salvando arquivo
end

clc
clear all
close all

lena = imread('lenaTest2.jpg');  %adquire imagem Lena como matriz
tam = size(lena); %obter tamanho da matriz da imagem
lena_v = double(reshape(lena,tam(1)*tam(2),1)); %reordena a matriz como um vetor
lena_b = de2bi(lena_v,8)'; %converte de decimal para binário
tam_b = size(lena_b); %obtem tamanho da nova matriz com dados em binário
lena_b = reshape(lena_b,[tam_b(1)*tam_b(2) 1]); %ajusta os bits para ficarem todos em uma coluna
tam_b = size(lena_b); %adquire as novas dimensões do vetor binário
pe = .25; %probabilidade de erro
ber = []; %pré alocando vetor com resultados ber
for n = 2:2:16 %laço para simulaçoes aumentando gradualmente o número de bits repetidos
    lena_n = repmat(lena_b,[1 n+1]); %repete os bits
    tam_n = size(lena_n); %adquire novas dimensões
    bsc_im = bsc(lena_n,pe); %simula canal bsc para um erro pe
    v_erro = sum(bsc_im,2); %soma cada bloco de bits repetidos
    resp = lena_n(:,1); %pré alocando resposta
    for p = 1:tam_n(1) %laço para corrigir os bits
        if v_erro(p) > (n+1)/2 %se a maioria dos bits de cada bloco for 1, resposta igual a 1
            resp(p,1) = 1;
        else %se a maioria dos bits de cada bloco for 0, resposta igual a 0
            resp(p,1) = 0;
        end
    end
    ber = [ber sum(abs(lena_b-resp),1)/tam_b(1)*100]; %calcula a ber
    figure %plotando imagem resultante
    im = uint8(bi2de(reshape(resp,[8 (tam_n(1))/8])')); %colocando no formato original
    im = reshape(im,tam); %colocando no formato original
    imshow(im) 
    xlabel(['n = ' num2str(n) ', ' '  BER = ' num2str(ber(end)) '$\%$' ' e $p_e=25\%$' ],'Interpreter','Latex','FontSize',14) %legenda
    saveas(gcf,['pe25' 'n' num2str(n) 'lena.pdf']) %salvando arquivo
end
figure %plotando ber
loglog(2:2:16,ber./100)
ylabel('BER','Interpreter','Latex','FontSize',14)
xlabel('$n$','Interpreter','Latex','FontSize',14)
grid on
saveas(gcf,'berxvsn.pdf')

clc
clear all
close all

lena = imread('lenaTest2.jpg'); %adquire imagem Lena como matriz
tam = size(lena); %obter tamanho da matriz da imagem
lena_v = double(reshape(lena,tam(1)*tam(2),1)); %reordena a matriz como um vetor
lena_b = de2bi(lena_v,8)';  %converte de decimal para binário
tam_b = size(lena_b); %obtem tamanho da nova matriz com dados em binário
lena_b = reshape(lena_b,[tam_b(1)*tam_b(2) 1]); %ajusta os bits para ficarem todos em uma coluna
tam_b = size(lena_b); %adquire as novas dimensões do vetor binário
pe = .25;  %probabilidade de erro
C = 1 + pe*log(pe)/log(2) + (1-pe)*log(1-pe)/log(2); %calcual a capacidade do canal BSC
n = ceil(1/C); %determina n para R_c = 1/n
tax = 1/n; %determina R_c máxima código ideal
lena_n = repmat(lena_b,[1 n+1]); %repetindo bits
tam_n = size(lena_n); %adquire novo tamanho dos dados
bsc_im = bsc(lena_n,pe); %obtem resultado após passar pelo canal BSC com erro pe = 25%
v_erro = sum(bsc_im,2);  %soma cada bloco de bits repetidos
resp = lena_n(:,1); %pré alocando resposta
for p = 1:tam_n(1) %laço para corrigir os bits
    if v_erro(p) > (n+1)/2 %se a maioria dos bits de cada bloco for 1, resposta igual a 1
        resp(p,1) = 1;
    else %se a maioria dos bits de cada bloco for 0, resposta igual a 0
        resp(p,1) = 0;
    end
end
ber = sum(abs(lena_b-resp),1)/tam_b(1)*100; %calcula ber  
figure %gera imagem
im = uint8(bi2de(reshape(resp,[8 (tam_n(1))/8])')); %colaca imagem no formato original 
im = reshape(im,tam); %colaca imagem no formato original 
imshow(im)
xlabel(['n = ' num2str(n) ', ' '  BER = ' num2str(ber(end)) '$\%$' ' e $p_e = 25\%$'],'Interpreter','Latex','FontSize',18) %legenda
saveas(gcf,['nidealpe25lena.pdf']) %salvando arquivo
ber_op = 999; %pré alocando ber para verificação manual
n = 2; %realocando valor de repetições
while ber_op > 10^(-5) %laço para determinar solução com R_c = 1/n com BER < 10^(-5)
    lena_n = repmat(lena_b,[1 n+1]); %repetindo bits
    tam_n = size(lena_n); %adquire novo tamanho dos dados
    bsc_im = bsc(lena_n,pe); %obtem resultado após passar pelo canal BSC com erro pe = 25%
    v_erro = sum(bsc_im,2);  %soma cada bloco de bits repetidos
    resp = lena_n(:,1); %pré alocando resposta
        for p = 1:tam_n(1) %laço para corrigir os bits
            if v_erro(p) > (n+1)/2 %se a maioria dos bits de cada bloco for 1, resposta igual a 1
                resp(p,1) = 1;
            else %se a maioria dos bits de cada bloco for 0, resposta igual a 0
                resp(p,1) = 0;
            end
        end
    ber_op = sum(abs(lena_b-resp),1)/tam_b(1) %calcula ber
    n = n + 2 %incrementa número de repetições
end
figure %gera imagem
im = uint8(bi2de(reshape(resp,[8 (tam_n(1))/8])')); %colaca imagem no formato original 
im = reshape(im,tam); %colaca imagem no formato original 
imshow(im)
xlabel(['n = ' num2str(n-2) ', ' '  BER = ' num2str(ber_op*100) '$\%$' ' e $p_e = 25\%$'],'Interpreter','Latex','FontSize',18) %legenda
saveas(gcf,['n' num2str(n-2) 'pe25lena.pdf']) %salvando arquivo