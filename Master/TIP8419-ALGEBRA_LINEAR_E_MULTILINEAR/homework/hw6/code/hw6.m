%% [TIP8419 - Algebra Linear e Multilinear] - Homework 6
% by Lucas Abdalah
% ----------------------------
% 
% hw6.m
% Author: Lucas Abdalah
%
clearvars; close all;
clc; pause(0.1);


%% Unfold Example ----------------------------
% Example 2.6 Multi-way Analysis With Applications in the Chemical Sciences (Smilde, 2004)
I = 4; J = 3; K = 2;

X(:, :, 1) = reshape(([1:9, flip(1:3)].'), [J I]).';
X(:, :, 2) = reshape(([5:9, 4:5, 3, 2, 4:6].'), [J I]).';
X_Nmode = {nd.unfold(X, 1), nd.unfold(X, 2), nd.unfold(X, 3)};

file = 'hw6_unfold.txt';
nd.tensor2txt(file, X, 'w', 'Tensor X');
nd.mat2txt(file, X_Nmode{1}, 'a', 'Tensor X (mode-1)');
nd.mat2txt(file, X_Nmode{2}, 'a', 'Tensor X (mode-2)');
nd.mat2txt(file, X_Nmode{3}, 'a', 'Tensor X (mode-3)');


%% Fold Example ----------------------------
% Example 2.6 Multi-way Analysis With Applications in the Chemical Sciences (Smilde, 2004)
X_build1 = nd.fold(X_Nmode{1}, [I J K], 1);
X_build2 = nd.fold(X_Nmode{2}, [I J K], 2);
X_build3 = nd.fold(X_Nmode{3}, [I J K], 3);

file = 'hw6_fold.txt';
nd.mat2txt(file, X_Nmode{1}, 'w', 'Tensor X (mode-1)');
nd.mat2txt(file, X_Nmode{2}, 'a', 'Tensor X (mode-2)');
nd.mat2txt(file, X_Nmode{3}, 'a', 'Tensor X (mode-3)');
nd.tensor2txt(file, X_build1, 'a', 'Tensor X from (mode-1)');
nd.tensor2txt(file, X_build2, 'a', 'Tensor X from (mode-2)');
nd.tensor2txt(file, X_build3, 'a', 'Tensor X from (mode-3)');


%% Unfold Validation ----------------------------
d = load('Practice_5_unfolding_folding.mat');
d.X_Nmode = {nd.unfold(d.tenX, 1), nd.unfold(d.tenX, 2), nd.unfold(d.tenX, 3)};

fprintf('\n--------------- Unfold Validation -------------- \n')
fprintf('\tsum(X1 - unfold(X, 1)) = %2.2f \n\tsum(X2 - unfold(X, 2)) = %2.2f \n\tsum(X3 - unfold(X, 3)) = %2.2f \n', ...
    sum(d.X1 - d.X_Nmode{1}, 'all'), sum(d.X2 - d.X_Nmode{2}, 'all'), sum(d.X3 - d.X_Nmode{3}, 'all'))


%% Fold Example ----------------------------
[d.I, d.J, d.K] = size(d.tenX);

d.X_build1 = nd.fold(d.X1, [d.I, d.J, d.K], 1);
d.X_build2 = nd.fold(d.X2, [d.I, d.J, d.K], 2);
d.X_build3 = nd.fold(d.X3, [d.I, d.J, d.K], 3);

fprintf('\n--------------- Fold Validation -------------- \n')
fprintf('\tsum(tenX - fold(X1)) = %2.2f \n\tsum(tenX - fold(X2)) = %2.2f \n\tsum(tenX - fold(X3)) = %2.2f \n', ...
    sum(d.tenX - d.X_build1, 'all'), sum(d.tenX - d.X_build3, 'all'), sum(d.tenX - d.X_build3, 'all'))


%% N-mode Validation ----------------------------
d_ = load('Practice_5_product_mode_n.mat');
d_.tenY_test = nd.N_mode(d_.tenX,{d_.Z}, 1);
[~, Y_NMSEdB1] = nd.nmse(d_.tenY, d_.tenY_test);

fprintf('\n--------------- N-mode Product Validation -------------- \n')
fprintf('\tNMSE between a given tensor and its version afected by the 1-mode product: %2.2f dB\n', Y_NMSEdB1);