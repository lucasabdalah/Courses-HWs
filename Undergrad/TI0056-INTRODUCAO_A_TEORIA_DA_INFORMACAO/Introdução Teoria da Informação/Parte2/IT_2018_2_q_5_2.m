lena = imread('lenaTest2.jpg');
tam = size(lena);
lena_v = double(reshape(lena,tam(1)*tam(2),1));
lena_b = de2bi(lena_v,8)';
tam_b = size(lena_b);
lena_b = reshape(lena_b,[tam_b(1)*tam_b(2) 1]);
tam_b = size(lena_b);
pe = .25;
ber = [];
for n = 2:2:16
    lena_n = repmat(lena_b,[1 n+1]);
    tam_n = size(lena_n);
    bsc_im = bsc(lena_n,pe);
    v_erro = sum(bsc_im,2);
    resp = lena_n(:,1);
    for p = 1:tam_n(1)
        if v_erro(p) > (n+1)/2
            resp(p,1) = 1;
        else
            resp(p,1) = 0;
        end
    end
    ber = [ber sum(abs(lena_b-resp),1)/tam_b(1)*100];
    figure
    im = uint8(bi2de(reshape(resp,[8 (tam_n(1))/8])'));
    im = reshape(im,tam);
    imshow(im)
    xlabel(['n = ' num2str(n) ', ' '  BER = ' num2str(ber(end)) '$\%$' ' e $p_e=25\%$' ],'Interpreter','Latex','FontSize',14)
    saveas(gcf,['n' num2str(n) 'pe25lena.pdf'])
end
figure
loglog(2:2:16,ber./100)
ylabel('BER','Interpreter','Latex','FontSize',14)
xlabel('$n$','Interpreter','Latex','FontSize',14)
grid on
saveas(gcf,'ber_x_n.pdf')
