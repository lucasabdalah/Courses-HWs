lena = imread('lenaTest2.jpg');
tam = size(lena);
lena_v = double(reshape(lena,tam(1)*tam(2),1));
lena_b = de2bi(lena_v,8);
tam_b = size(lena_b);
pe = 0:.25:1;
bsc_im = zeros([tam_b length(pe)]);
for p = 1:length(pe)
    bsc_im(:,:,p) = bsc(lena_b,pe(p));
end
for p = 1:length(pe)
    im = uint8(reshape(bi2de(bsc_im(:,:,p)),tam));
    figure
    imshow(im)
    xlabel(['$p_e = ' , num2str(pe(p)*100), ' \%$'],'Interpreter','Latex','FontSize',20)
    saveas(gcf,[num2str(pe(p)*100) 'errolena.pdf'])
end