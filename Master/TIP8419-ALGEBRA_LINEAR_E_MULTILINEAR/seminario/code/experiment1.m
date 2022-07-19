%% [TIP8419 - Algebra Linear e Multilinear] - Seminar
% Author: Lucas Abdalah
% ----------------------------
% Partial reproduction of Figure 1: 
%   Alternating group lasso for block-term tensor decomposition with application 
%   to ECG source separation 
%   J. H. de M. Goulart et al., IEEE TSP, 2020
%
% experiment1.m
%

clearvars; close all;


%% General Setup
Ls = [6 5 4];
Is = [18 18 4];
R = length(Ls);
Le = 6;
min_dist = 0.1;
SNR = 20;
RMC = 500;


%% Tensorlab - Parameters
opts_tlab.Display = false;
opts_tlab.Compression = false;
opts_tlab.TolX = 1e-5;
opts_tlab.TolF = 1e-5;
opts_tlab.MaxIter = 1.5e3;


%% Tensorlab True Ranks - Save vector
err_Y_TL_TR = zeros(RMC, 1);
err_blocks_TL_TR = zeros(RMC, 3);
R_TL_TR = zeros(RMC, 3);
t_TL_TR = zeros(RMC, 1);
it_TL_TR = zeros(RMC, 1);
Xcols_TL_TR = cell(RMC, 1);
Perms_TL_TR = zeros(RMC, 3);
Rnd_init_TL_TR = zeros(RMC, 1);


%% Tensorlab Same Ranks - Save vector
err_Y_TL_SR = zeros(RMC, 1);
err_blocks_TL_SR = zeros(RMC, 3);
R_TL_SR = zeros(RMC, 3);
t_TL_SR = zeros(RMC, 1);
it_TL_SR = zeros(RMC, 1);
Xcols_TL_SR = cell(RMC, 1);
Rnd_init_TL_SR = zeros(RMC, 1);


%% Monte Carlo simulation

filename = '';

for rmc = 1:RMC

    %% Generating model
    A = randn(Is(1), sum(Ls));
    B = randn(Is(2), sum(Ls));
    X = normcol(randn(Is(3), R));
    while min(vec((1 - abs(X' * X) + diag(inf * ones(1,R))))) <= min_dist
        X = normcol(randn(Is(3), R));
    end

    [Y0, Yb] = ll1btd(A, B, X, Ls);

        %% Initial point
        A0 = randn(Is(1), R*Le);
        B0 = randn(Is(2), R*Le);
        X0 = randn(Is(3), R);
        X0 = normcol(X0);
        Y = noisybtd(Y0, SNR); % Generating noisy data

        %% Tensorlab with same ranks --------------------------------------

        disp(' ')

        exttimer = tic;
        [U, output] = ll1(Y, Le*[1 1 1], opts_tlab);
        t_TL_SR(rmc) = toc(exttimer);
        it_TL_SR(rmc) = output.Algorithm.iterations;

        Rnd_init_TL_SR(rmc) = output.Initialization.Backup;

        Ae = zeros(Is(1), R*Le);
        Be = zeros(Is(2), R*Le);
        Xe = zeros(Is(3), R);
        
        for r = 1:R
            inds = (r-1)*Le+1 : r*Le;
            Ae(:,inds) = U{r}{1};
            Be(:,inds) = U{r}{2} * U{r}{4};
            Xe(:,r) = U{r}{3};
            R_TL_SR(rmc, r) = rank(Ae(:,inds) * Be(:,inds).');
        end

        [Ye, Ybe] = ll1btd(Ae, Be, Xe, Le * ones(1,R));
        [ps, err_blocks_TL_SR(rmc, :)] = permblock(Yb, Ybe);
        fprintf('Error Tensorlab - same ranks: %4.2e  %4.2e  %4.2e, mean = %4.2e     (%5d iterations) \n',...
            err_blocks_TL_SR(rmc, :),mean(err_blocks_TL_SR(rmc, :)),it_TL_SR(rmc));
        err_Y_TL_SR(rmc) = norm(Y0(:) - Ye(:))^2/norm(Y0(:))^2;

        Xcols_TL_SR{rmc} = normcol(Xe);

        %% Tensorlab with true ranks --------------------------------------

        exttimer = tic;
        [U, output] = ll1(Y, Ls, opts_tlab);
        t_TL_TR(rmc) = toc(exttimer);
        it_TL_TR(rmc) = output.Algorithm.iterations;

        Rnd_init_TL_TR(rmc) = output.Initialization.Backup;

        Ae = zeros(Is(1), sum(Ls));
        Be = zeros(Is(2), sum(Ls));
        Xe = zeros(Is(3), R);
        for r = 1:R
            inds = sum(Ls(1:r-1))+1:sum(Ls(1:r));
            Ae(:,inds) = U{r}{1};
            Be(:,inds) = U{r}{2} * U{r}{4};
            Xe(:,r) = U{r}{3};
            R_TL_TR(rmc, r) = rank(Ae(:,inds) * Be(:,inds).');
        end

        [Ye, Ybe] = ll1btd(Ae, Be, Xe, Ls);
        [ps, err_blocks_TL_TR(rmc, :)] = permblock(Yb, Ybe);
        fprintf('Error Tensorlab - true ranks: %4.2e  %4.2e  %4.2e, mean = %4.2e     (%5d iterations) \n',...
            err_blocks_TL_TR(rmc, :),mean(err_blocks_TL_TR(rmc, :)),it_TL_TR(rmc));
        err_Y_TL_TR(rmc) = norm(Y0(:) - Ye(:))^2/norm(Y0(:))^2;
        Perms_TL_TR(rmc,:) = ps;
        disp(['Optimal permutation = ' num2str(ps)])

        Xcols_TL_TR{rmc} = normcol(Xe);

    try
        delete([filename '.mat']);
    catch
    
    end
    filename = 'MC_scen1_rmc500';
    save([filename '.mat']);
    % datetimestamp = datestr(now, 'yyyy-mm-dd');
    % save([filename '.mat']);

    disp(['Realização ' num2str(rmc) ' concluída. ============================= '])
    disp(' ')

end

%% Plot results

figure
ecdf(vec(err_blocks_TL_SR))
hold on
ecdf(vec(err_blocks_TL_TR))
xlim([0 1])
set(gca,'xscale','log')
legend('Tensorlab (same ranks)','Tensorlab (true ranks)')
grid on

sum((Perms_TL_TR(:,1) == 1) & (Perms_TL_TR(:,2) == 2) & (Perms_TL_TR(:,3) == 3))/RMC

sum(Rnd_init_TL_SR)/RMC

sum(Rnd_init_TL_TR)/RMC

%% Local Functions ------------------------------------------------------------

function [Xn,ns] = normcol(Xn)
    
    ns = sqrt(sum(abs(Xn).^2));
    
    i0 = (ns == 0);
    
    ns(i0) = 1;
    
    Xn = Xn ./ repmat( ns, [size(Xn,1) 1] );
    
    ns(i0) = 0;
    
end


function v = vec(V)
    v = V(:);
end


function [Yn, Nois, sigY] = noisybtd(Y, SNR)

    Nois = randn(size(Y));
    
    sigY = 10^(-SNR/20) * norm(Y(:))/norm(Nois(:));
    
    Yn = Y + sigY * Nois;
end