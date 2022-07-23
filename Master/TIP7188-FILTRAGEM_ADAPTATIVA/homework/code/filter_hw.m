%% MÉTODOS 
% [TIP7188 - Filtragem Adaptativa]
% Author: Lucas Abdalah
%
% filter_hw.m
% 
% filter_hw is a package developped for the Adaptative Filtering Course
% It is a way to make a compilation for all function
%
% CONTENT
%   
%   HOMEWORK 2 - PROBLEM 5
% 
%   SAVE DATA TO TXT FILE 
%       filter_hw.MAT2TXT      - Write a matrix X into a txt file
%       filter_hw.TENSOR2TXT   - Write a 3D tensor X into a txt file
% 
%   PLACE HOLDER
% 

classdef filter_hw

methods(Static)


%% HOMEWORK 2 - PROBLEM 5
function hw2p5(varargin)
% FILTER_HW.HW2P5  Perfom the error surface propose on the Hw 2, problem 5
%
%   See also.

    if isempty(varargin)
        save_results = false;
    else
        save_results = varargin{1};
    end

    N = 25;
    w_lim = 100;
    w = [linspace(-w_lim,w_lim,N); linspace(-w_lim,w_lim,N)];
    [w_0, w_1] = meshgrid(w(1,:), w(2,:));
    J_surface = @(w_0, w_1) 24.40 - 4.*w_0 - 9.*w_1 + w_0.^2 + w_1.^2;
    J = J_surface(w_0, w_1);
    h = figure();
    surf(w_0, w_1, J, 'EdgeColor', 'none');
    colormap turbo;
    xlabel('$w_0$', 'FontSize', 16, 'interpreter', 'latex');
    ylabel('$w_1$', 'FontSize', 16, 'interpreter', 'latex');
    zlabel('$J$', 'FontSize', 16, 'interpreter', 'latex');
    view([-24.5036297640653 47.6514617014408]);
    colorbar('box', 'off');
    grid on;
    axis tight;
    pathName = 'figures/';
    filter_hw.export_fig(save_results, h, [pathName, 'hw2p5']);
end


%% HOMEWORK 3 - PROBLEM 4

function filter_path(signal_d_var, weights, wiener, Rx, p, c_)
% FILTER_HW.FILTER_PATH  Perfom the weights path surface
%
%   See also.
    step = 0.25;
    X = meshgrid (-1:step:1,-1:step:1);
    w = [X(:), reshape(transpose(X),[],1)];
    [wLen, ~] = size(w);
    J = zeros(wLen, 1);
    for n = 1:wLen
        J(n) = signal_d_var - 2*w(n,:)*p + w(n,:)*Rx*w(n,:).';
    end
    contour(X, X', reshape(J,size(X)), '-.', 'color', 'k');
    hold on;
    scatter(weights(1,:), weights(2, :), '.', 'MarkerEdgeColor', c_);
    hold off;   
    ha = annotation('textarrow', [0 0], [0 0],'String', 'Wiener');
    ha.Parent = gca;
    ha.X = [wiener(1)+0.15 wiener(1)];
    ha.Y = [wiener(2)-0.4 wiener(2)];
    grid on;
    
end

function [error, weights] = dga(signal_x, signal_d, order, mi, Rx, p)
% FILTER_HW.DGA  Perfom the Deterministic Gradient Algorithm 
%
%   See also.
    N = length(signal_x);
    error = zeros(N,1);
    weights = zeros(order, N);
    signal_d = signal_d(order:end,1);
    
    for n = 1:(N - order - 1)
        error(n,1) = signal_d(n) - weights(:,n)'*signal_x(n:n+order-1);
        weights(:,n+1) = weights(:,n) - 2*mi*(Rx*weights(:,n) - p);
    end
end


function [error, weights] = lms(signal_x, signal_d, order, mi)
% FILTER_HW.LMS  Perfom the LMS Algorithm 
%
%   See also.
    N = length(signal_x);
    error = zeros(N,1);
    weights = zeros(order, N);
    signal_d = signal_d(order:end,1); 

    for n = 1:(N - order - 1)
        error(n) = signal_d(n) - weights(:,n)' * signal_x(n:n+order-1); 
        weights(:,n+1) = weights(:,n) + 2 * mi * error(n) * signal_x(n:n+order-1);
    end
    weights = flip(weights);
end


function [error, weights] = newton(signal_x, signal_d, order, mi, wiener)
% FILTER_HW.NEWTON  Perfom the Newton Algorithm 
%
%   See also.
    N = length(signal_x);
    error = zeros(N,1);
    weights = zeros(order, N);
    signal_d = signal_d(order:end,1);

    for n = 1:(N - order - 1)
        error(n,1) = signal_d(n) - weights(:,n)'*signal_x(n:n+order-1);
        weights(:,n+1) = weights(:,n) - mi*(weights(:,n) - wiener);
    end
end


function [error, weights] = nlms(signal_x, signal_d, order, mi, gamma)
% FILTER_HW.NLMS  Perfom the NLMS Algorithm 
%
%   See also.
    N = length(signal_x);
    error = zeros(N,1);
    weights = zeros(order, N);
    signal_d = signal_d(order:end,1);

    for n = 1:(N - order - 1)
        mi_normalized = mi/(gamma + norm(signal_x));
        error(n) = signal_d(n) - weights(:,n)' * signal_x(n:n + order-1); 
        weights(:,n+1) = weights(:,n) +  2 * mi_normalized * error(n) * signal_x(n:n+order-1);
    end
    weights = flip(weights);
end


function hw3p4(varargin)

    % Save or not the results 
    if isempty(varargin)
        save_results = false;
    else
        save_results = varargin{1};
    end

    pathName = 'figures/';

    h0 = figure();
    viscircles([0, 0], 1,'Color','k', 'LineStyle','-', 'LineWidth', 1.5);
    line([0 0],[-1 1], 'Color', 'k', 'HandleVisibility','off');
    line([-1 1], [0 0], 'Color', 'k', 'HandleVisibility','off');
    hold on
    scatter(-1.6, 0, 'o', 'filled');
    scatter(0.45, 0, 'o', 'filled');
    hold off
    xlabel('$\Re (Z)$', 'interpreter', 'latex');
    ylabel('$\Im (Z)$', 'interpreter', 'latex');
    axis([-1.7 1.7 -1.7 1.7]);
    legend('Channel Zeros', 'Filter Zeros', 'Location', 'Northeast');
    grid minor
    axis square

    filter_hw.export_fig(save_results, h0, [pathName, 'hw3p4-zeros']);

    
    % Color scheme to plot ------------------------------------
    c_ = struct('dg', [57 106 177]./255, 'lms', [204 37 41]./255, 'newton', [62 150 81]./255, 'nlms', [107 76 154]./255, 'mean', 'k');
    % General Setup ------------------------------------
    N = 1000; % Number of samples
    order = 2; % Filter order
    % Signal Model ------------------------------------
    signal_d = randn(N,1);
    signal_d_var = var(signal_d);
    % Noisy Version ------------------------------------
    Hz = [1 1.6]; 
    signal_x = filter(Hz,1,signal_d); 
    noise = sqrt(1/(10^(inf/10))).*randn(N,1);
    signal_x = signal_x + noise; 
    % Wiener Filter ------------------------------------
    Rxcorr = sort(xcorr(Hz));
    Rx = reshape([Rxcorr(end) Rxcorr], [2, 2]); % Autocorrelation matrix
    p = eye(2,1); % Cross-correlation
    wiener = Rx\p; % Wiener solution
    fprintf('Wiener solution: %2.2f \n %2.2f \n', wiener);
    % Deterministic Gradient Algorithm ------------------------------------
    dg.mi = 1e-2;
    [dg.error, dg.weights] = filter_hw.dga(signal_x, signal_d, order, dg.mi, Rx, p);
    % Newton Implementation ------------------------------------
    newton.mi = 5e-2;
    [newton.error, newton.weights] = filter_hw.newton(signal_x, signal_d, order, newton.mi, wiener);
    % LMS Algorithm ------------------------------------
    lms.mi = 1e-3;
    [lms.error, lms.weights] = filter_hw.lms(signal_x, signal_d, order, lms.mi);
    % NLMS Algorithm ------------------------------------
    nlms.mi = 5e-1;
    gamma = 0.5;
    [nlms.error, nlms.weights] = filter_hw.nlms(signal_x, signal_d, order, nlms.mi, gamma);
    % Plot - Deterministic Gradient Algorithm ------------------------------------
    h1 = figure(1)
    subplot(2,1,1);
    semilogy(1:N, dg.error.^2,'-','color', c_.dg , "linewidth", 1); % MSE
    hold on
    semilogy(1:N, repelem(mean(dg.error.^2), N), '--', 'color', c_.mean, "linewidth", 1);
    hold off
    xlabel('Iterations');
    ylabel('MSE');
    legend('Deterministic Gradient', 'Mean', 'Location', 'Best')
    grid on;
    axis tight
    subplot(2,1,2);
    filter_hw.filter_path(signal_d_var, dg.weights, wiener, Rx, p, c_.dg); % Solution Path
    xlabel('$w_1$', 'interpreter', 'latex');
    ylabel('$w_0$', 'interpreter', 'latex');
    legend('Solution Contour', 'Deterministic Gradient', 'Location', 'Northeast')
    axis tight
    
    filter_hw.export_fig(save_results, h1, [pathName, 'hw3p4-dga']);

    % Plot - Newton Implementation ------------------------------------
    h2 = figure(2)
    subplot(2,1,1);
    semilogy(1:N, newton.error.^2,'-','color', c_.newton, "linewidth", 1); % MSE Curve
    hold on
    semilogy(1:N, repelem(mean(newton.error.^2), N), '--', 'color', c_.mean, "linewidth", 1);
    hold off
    xlabel('Iterations');
    ylabel('MSE');
    legend('Newton', 'Mean', 'Location', 'Best');
    grid on;
    axis tight
    subplot(2,1,2);
    filter_hw.filter_path(signal_d_var, newton.weights, wiener, Rx, p, c_.newton);
    xlabel('$w_1$', 'interpreter', 'latex');
    ylabel('$w_0$', 'interpreter', 'latex');
    legend('Solution Contour', 'Newton',  'Location', 'Northeast')
    axis tight
    
    filter_hw.export_fig(save_results, h2, [pathName, 'hw3p4-newton']);
    
    % Plot - LMS Algorithm ------------------------------------
    h3 = figure(3)
    subplot(2,1,1);
    semilogy(1:N, lms.error.^2,'-','color', c_.lms , "linewidth", 1); % MSE
    hold on
    semilogy(1:N, repelem(mean(lms.error.^2), N), '--', 'color', c_.mean, "linewidth", 1);
    hold off
    xlabel('Samples');
    ylabel('MSE');
    legend('LMS', 'Mean', 'Location', 'Best')
    grid on;
    axis tight
    subplot(2,1,2);
    filter_hw.filter_path(signal_d_var, lms.weights, wiener, Rx, p, c_.lms); % Solution Path
    xlabel('$w_1$', 'interpreter', 'latex');
    ylabel('$w_0$', 'interpreter', 'latex');
    legend('Solution Contour', 'LMS',  'Location', 'Northeast')
    axis tight
    
    filter_hw.export_fig(save_results, h3, [pathName, 'hw3p4-lms']);

    % Plot - NLMS Implementation ------------------------------------
    h4 = figure(4)
    subplot(2,1,1);
    semilogy(1:N, nlms.error.^2,'-','color', c_.nlms, "linewidth", 1);
    hold on
    semilogy(1:N, repelem(mean(nlms.error.^2), N), '--', 'color', c_.mean, "linewidth", 1);
    hold off
    xlabel('Samples');
    ylabel('MSE');
    legend('NLMS','Mean', 'Location', 'Best');
    grid on;
    axis tight
    subplot(2,1,2);
    filter_hw.filter_path(signal_d_var, nlms.weights, wiener, Rx, p, c_.nlms);
    xlabel('$w_1$', 'interpreter', 'latex');
    ylabel('$w_0$', 'interpreter', 'latex');
    legend('Solution Contour', 'NLMS',  'Location', 'Northeast')
    axis tight

    filter_hw.export_fig(save_results, h4, [pathName, 'hw3p4-nlms']);

end


%% ----------------------- OK ATE AQUI

%% HOMEWORK 3 - PROBLEM 5

function hw3p5(varargin)
% hw3p5
% Learning rate
    
    mi = 1/(97*2);
    % Filter order
    % I first implemented thinking of python notation, later I found out that
    % the reference book defines the order a bit different from what I usually
    % work. So to make the code close to Diniz notation the 'order + 1' is
    % needed.
    order = 15 + 1;
    % Number of samples
    Samples = 5000 + order;
    % Defining the mse error and filter coeficients vectors.
    error = zeros(Samples,1);
    weights = zeros(order, Samples);

    % Defining the energy of the noise vector as 1e-3.
    SNR_dB = 30;
    SNR_li = 10^(SNR_dB/10);
    variance_noise = 1/SNR_li;
    noise = sqrt(variance_noise).*randn(Samples,1);

    % Generating the original signal.
    signal_d = randn(Samples,1);
    signal_d = (signal_d-mean(signal_d))/std(signal_d);

    % Convolving the channel and the signal. To prevent the missmatch between 
    % the filtered signal and the desired signal we use filtfilt instead of filter or conv.
    Hz = [1 1 1 1 1 1 1 1 1 1 1 1];
    signal_x = filtfilt(Hz,1,signal_d);
    % Generating the noisy received signal.
    signal_x = signal_x + noise;
    signal_x = (signal_x-mean(signal_x))/std(signal_x);
    signal_d_hat = zeros(size(signal_d));
    for ss = 1:(Samples - order)
        signal_d_hat(ss) = weights(:,ss)'*signal_x(ss:ss+order-1);
        % Error between the desired signal and the filtered signal.
        error(ss) = signal_d(ss) - weights(:,ss)' * signal_x(ss:ss+order-1); 
        % Recursive expression.
        weights(:,ss+1) = weights(:,ss) + 2 * mi * error(ss) * signal_x(ss:ss+order-1);
    end
    signal_d_hat = (signal_d_hat-mean(signal_d_hat))/std(signal_d_hat);

    % MSE Curve
    figure
    semilogy(1:Samples, error.^2,'-','color', [0.3010 0.7450 0.9330], "linewidth", 1, "markersize", 8);
    title('LMS Behavior');
    xlabel('Samples');
    ylabel('MSE');
    grid on;
    % saveas(gcf,'L3Q5_mu_2.png')

    % Filter Response
    %https://www.mathworks.com/help/signal/ug/frequency-response.html#:~:text=To%20convert%20normalized%20frequency%20back,by%20half%20the%20sample%20frequency.&text=freqz%20can%20also%20accept%20a,(b%2Ca%2Cw)%3B
    [Hf,wf] = freqz(weights(:,ss + 1).',1,'whole',512);
    [Hc,wc] = freqz([1 1 1 1 1 1 1 1 1 1 1 1],1,'whole',512);
    txt = ['Unknown System'];
    plot(wc/pi,20*log10(abs(Hc)),'-','color', [0.3010 0.7450 0.9330], "linewidth", 1, "markersize", 8, "DisplayName", txt);
    hold on;
    txt = ['Filter Response'];
    plot(wf/pi,20*log10(abs(Hf)),'-','color', [0.4660 0.6740 0.1880], "linewidth", 1, "markersize", 8, "DisplayName", txt);
    title('System Identification with LMS')
    xlabel('Normalized Frequency (\times\pi rad/sample)')
    ylabel('Magnitude (dB)')
    grid on;
    legend_copy = legend("location", "southwest");
    set (legend_copy, "fontsize", 6);
    % saveas(gcf,'L3Q5_filter_response_2.png')

    % Temporal Evolution
    figure
    txt = ['Original Signal'];
    plot(1:Samples,signal_d,'-','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
    hold on;
    txt = ['Estimated Signal'];
    plot(1:Samples,signal_d_hat,'-','color', [0.4660 0.6740 0.1880], "linewidth", 2, "markersize", 8, "DisplayName", txt);
    title('Temporal Evolution');
    xlabel('Samples');
    xlim([0 250]);
    ylabel('Magnitude');
    grid on;
    legend_copy = legend("location", "southwest");
    set (legend_copy, "fontsize", 6);
    % saveas(gcf,'L3Q5_t.png')

    % Filter Response (Another method)
    % https://www.mathworks.com/matlabcentral/answers/514720-how-to-use-freqz-to-plot-filter-frequency-response
    %figure
    %Samples = 5e3;
    %fhz = 1e8;
    %[h_filter, w_filter] = freqz(weights(:,ss+1).',1,Samples,fhz);
    %[h_channel, w_channel] = freqz([1 1 1 1 1 1 1 1 1 1 1 1],1,Samples,fhz);
    %[h_channel, w_channel] = freqz([1 1.6],[1],Samples,fhz);
    %txt = ['Filter Response'];
    %plot(w_filter,abs(h_filter),'-','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
    %hold on;
    %txt = ['System Response'];
    %plot(w_channel,abs(h_channel),'-','color', [0.4660 0.6740 0.1880], "linewidth", 2, "markersize", 8, "DisplayName", txt);
    %hold off;
    %legend_copy = legend("location", "southwest");
    %set (legend_copy, "fontsize", 12);
    %grid on;
    %set(gca,'XScale','log')
    %set(gca,'YScale','log')
    %title('Filter Response')
    %xlabel('Frequency (Hz)')
    %ylabel('Filter Response')
    % saveas(gcf,'L3Q5_filter_response.png')


    % The filter order is 12
    order = 12;
    Rxx = zeros(order,order);
    for k = 1:10000
        x = randn(10000,1) + randn(10000,1);
        % It is necessary to centralize the data.
        x = (x-mean(x))/std(x);
        y = zeros(length(x) + order - 1,1);
        for i = order:length(x)
        y(i + order - 1) = x(i) + x(i-1) + x(i-2) +  x(i-3) + x(i-4) + x(i-5) + x(i-6) +...
            x(i-7) + x(i-8) + x(i-9) + x(i-10) + x(i-11);
        %y(i-1) = x(i) + 1.6*x(i-1);
        end

        % The following references were used to solve this problem.
        %https://www.mathworks.com/matlabcentral/answers/4580-calculation-of-autocorrelation-matrix
        %https://ccrma.stanford.edu/~jos/sasp/Sample_Autocorrelation.html
        %http://matlab.izmiran.ru/help/toolbox/signal/corrmtx.html
        %https://www.mathworks.com/help/signal/ref/corrmtx.html
        [~,R] = corrmtx(y,order - 1,'autocorrelation');
        Rxx = Rxx + R;
    end
    Rxx = Rxx./10000;
    step  = 1/max(eig(Rxx)) 

    % Misadjustment for mu/2, mu/10 and mu/50
    Rxx_empirico = Rxx;

    %     misadjustment_e_2 = ((0.05/2)*(trace(Rxx_empirico)))/(1 - (0.05/2)*(trace(Rxx_empirico)));
%     misadjustment_e_10 = ((0.05/10)*(trace(Rxx_empirico)))/(1 - (0.05/10)*(trace(Rxx_empirico)));
%     misadjustment_e_50 = ((0.05/50)*(trace(Rxx_empirico)))/(1 - (0.05/50)*(trace(Rxx_empirico)));
%     Rxx_teorico = ceil(Rxx);
%     misadjustment_t_2 = ((0.05/2)*(trace(Rxx_teorico)))/(1 - (0.05/2)*(trace(Rxx_teorico)));
%     misadjustment_t_10 = ((0.05/10)*(trace(Rxx_teorico)))/(1 - (0.05/10)*(trace(Rxx_teorico)));
%     misadjustment_t_50 = ((0.05/50)*(trace(Rxx_teorico)))/(1 - (0.05/50)*(trace(Rxx_teorico)));
% end


% %% HOMEWORK 3 - PROBLEM 6

% function hw3p6(varargin)
%     % hw3p6
%     % (a) --------------------------------------
%     disp('a')
%     % Training Phase
%     % Learning rate
%     mi = 0.4;
%     gamma = 1e-3;
%     % I first implemented thinking of python notation, later I found out that
%     % the reference book defines the order a bit different from what I usually
%     % work. So to make the code close to Diniz notation the 'order + 1' is
%     % needed.
%     order = 15 + 1;
%     % Number of samples
%     Samples = 500;
%     % Defining the mse error and filter coeficients vectors.
%     error = zeros(Samples,1);
%     weights = zeros(order, Samples);

%     % Defining the energy of the noise vector.
%     SNR = inf;
%     QAM_train = 4;
%     signal_d_train = randi([0,QAM_train - 1],[Samples 1]); % The same pilot for every pilot frame and block.
%     signal_d_train = qammod(signal_d_train,QAM_train); % 4-QAM Pilot Signal.

%     % Convolving the channel and the signal.It is necesssary to sincronize the signal.
%     %https://www.mathworks.com/help/signal/ref/filtfilt.html
%     Hz = [0.5 1.2 1.5 -1];
%     signal_x_train = filtfilt(Hz,1,signal_d_train);

%     % Training noise
%     snr = 10^(SNR/10);
%     energy_symbol = mean(abs(signal_x_train(:)).^2); % Energy symbol pilot. 
%     var_noise = energy_symbol .*  1/snr; % Variance of the noise.
%     noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

%     % Generating the noisy received signal.
%     signal_x_train = signal_x_train + noise;
%     % NLMS algorithm
%     for s = order:Samples
%         aux = signal_x_train(s:-1:s-order+1);
%         mi_normalized = mi/(gamma + norm(aux)^2);
%         error(s) = signal_d_train(s-order+1) - weights(:,s)'*aux;
%         % Recursive expression.
%         weights(:,s+1) = weights(:,s) + mi_normalized * conj(error(s)) * aux;
%     end

%     % Transmission Phase
%     % Number of samples
%     Samples = 5000 + order;
%     % Defining the mse error and filter coeficients vectors.
%     aux = weights(:,s+1);
%     error = zeros(Samples,1);
%     weights = zeros(order, Samples);
%     weights(:,order) = aux;

%     % Defining the energy of the noise vector.
%     SNR = 30;
%     QAM = 16;
%     signal_d = randi([0,QAM - 1],[Samples 1]); % The same pilot for every pilot frame and block.
%     signal_d = qammod(signal_d,QAM); % 4-QAM Pilot Signal.

%     % Convolving the channel and the signal. 
%     signal_x = filtfilt(Hz,1,signal_d);

%     % Transmission noise
%     snr = 10^(SNR/10);
%     energy_symbol = mean(abs(signal_x(:)).^2); % Energy symbol pilot. 
%     var_noise = energy_symbol .*  1/snr; % Variance of the noise.
%     noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

%     % Generating the noisy received signal.
%     signal_x = signal_x + noise;
%     signal_d_hat = zeros(size(signal_d));
%     % NLMS algorithm
%     for s = order:Samples
%         aux = signal_x(s:-1:s-order+1);
%         mi_normalized = mi/(gamma + norm(aux)^2);
%         % Filtering the signal
%         signal_d_hat(s-order+1) = weights(:,s)'*aux;
%         % The equalizer does know the original signal
%         %error(s) = signal_d(s-order+1) - weights(:,s)'*aux;
%         % The equalizer does not know the original signal
%         error(s) = qammod(qamdemod(signal_x(s-order+1),QAM),QAM) - weights(:,s)'*aux;
%         % Recursive expression.
%         weights(:,s+1) = weights(:,s) + mi_normalized * conj(error(s)) * aux;
%     end

%     % MSE Curve
%     figure
%     semilogy(1:Samples, abs(error).^2,'-','color', [0.3010 0.7450 0.9330], "linewidth", 1, "markersize", 8);
%     title('NLMS Behavior');
%     xlabel('Samples');
%     xlim([0 5000]);
%     ylabel('MSE');
%     grid on;
%     % saveas(gcf,'L3_A_mse.png')

%     % Temporal Evolution
%     aux = qamdemod(signal_d_hat,QAM);
%     aux1 = aux(1:100);
%     aux2 = aux(4900:5000);

%     figure
%     subplot(211)
%     txt = ['Original Signal'];
%     stem(1:100, qamdemod(signal_d(1:100),QAM),'-','color', [0.3010 0.7450 0.9330], "linewidth", 1, "markersize", 3, "DisplayName", txt);
%     hold on;
%     txt = ['Estimated Signal'];
%     stem(1:100, aux1,'-','color', [0.4660 0.6740 0.1880], "linewidth", 1, "markersize", 3, "DisplayName", txt);
%     hold off;
%     title('Temporal Evolution of First Samples');
%     xlabel('Sample');
%     ylabel('Magnitude');
%     legend_copy = legend("location", "southwest");
%     set (legend_copy, "fontsize", 6);
%     grid on;
%     subplot(212)
%     txt = ['Original Signal'];
%     stem(4900:5000, qamdemod(signal_d(4900:5000),QAM),'-','color', [0.3010 0.7450 0.9330], "linewidth", 1, "markersize", 3, "DisplayName", txt);
%     hold on;
%     txt = ['Estimated Signal'];
%     stem(4900:5000, aux2,'-','color', [0.4660 0.6740 0.1880], "linewidth", 1, "markersize", 3, "DisplayName", txt);
%     hold off;
%     title('Temporal Evolution of Last Samples');
%     xlabel('Sample');
%     ylabel('Magnitude');
%     legend_copy = legend("location", "southwest");
%     set (legend_copy, "fontsize", 6);
%     grid on;
%     % saveas(gcf,'L3_A_t.png')

%     % Constellation
%     % https://www.mathworks.com/help/comm/gs/examine-16-qam-using-matlab.html
%     figure
%     subplot(221)
%     plot(signal_d_train,'.','color', 'y',"markersize", 8)
%     title('Training Signal');
%     xlabel('In Phase');
%     ylabel('Quadrature');
%     axis([-2 2 -2 2]);
%     set(gca,'Color','k');
%     grid on;
%     subplot(222)
%     plot(signal_d,'.','color', 'y',"markersize", 8)
%     title('Original Signal');
%     xlabel('In Phase');
%     ylabel('Quadrature');
%     set(gca,'Color','k');
%     grid on;
%     subplot(223)
%     plot(signal_x,'.','color', 'y',"markersize", 8)
%     title('Transmitted Signal');
%     xlabel('In Phase');
%     ylabel('Quadrature');
%     set(gca,'Color','k');
%     grid on;
%     subplot(224)
%     plot(qammod(qamdemod(signal_d_hat,QAM),QAM),'.','color', 'y',"markersize", 8)
%     title('Filtered Signal + Decisor');
%     xlabel('In Phase');
%     ylabel('Quadrature');
%     set(gca,'Color','k');
%     grid on;
%     % saveas(gcf,'L3_A_c.png')
    
%     % (b) --------------------------------------
%     disp('b')
%     % Learning rate
%     mi = 1e-3;
%     % Filter order
%     order = 15 + 1;

%     % Generating Transmisssion Phase Data
%     % If we want to evaluate the impact of the training size then we should
%     % apply the same transmitted sequence to all the cases
%     % Number of samples
%     Samples = 5000 + 50;
%     % Defining the energy of the noise vector.
%     SNR = 30;
%     QAM = 16;
%     signal_d = randi([0,QAM - 1],[Samples 1]); 
%     signal_d = qammod(signal_d,QAM);
%     % Convolving the channel and the signal.
%     Hz = [0.5 1.2 1.5 -1];
%     signal_x = filter(Hz,1,signal_d);
%     % Transmision noie
%     snr = 10^(SNR/10);
%     energy_symbol = mean(abs(signal_x(:)).^2); % Energy symbol pilot. 
%     var_noise = energy_symbol .*  1/snr; % Variance of the noise.
%     noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));
%     % Generating the noisy received signal.
%     signal_x = signal_x + noise;

%     % Training with 50 Samples

%     % Training Phase

%     % Number of samples
%     Samples = 50;

%     % Defining the mse error and filter coeficients vectors.
%     error = zeros(Samples,1);
%     weights = zeros(order, Samples);

%     % Defining the energy of the noise vector.
%     QAM_train = 4;
%     signal_d_train = randi([0,QAM_train - 1],[Samples 1]);
%     signal_d_train = (1/sqrt(2)) * qammod(signal_d_train,QAM_train); 

%     % Convolving the channel and the signal.
%     Hz = [0.5 1.2 1.5 -1];
%     signal_x_train = filter(Hz,1,signal_d_train);

%     % Training noise
%     snr = 10^(inf/10);
%     energy_symbol = mean(abs(signal_x_train(:)).^2); % Energy symbol pilot. 
%     var_noise = energy_symbol .*  1/snr; % Variance of the noise.
%     noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

%     % Generating the noisy received signal.
%     signal_x_train = signal_x_train + noise;
%     % LMS algorithm
%     for s = order:Samples
%         aux = signal_x_train(s:-1:s-order+1);
%         error(s) = signal_d_train(s-order+1) - weights(:,s)'*aux;
%         % Recursive expression.
%         weights(:,s+1) = weights(:,s) + 2 * mi * conj(error(s)) * aux;
%     end

%     % Transmission Phase

%     % Defining the mse error and filter coeficients vectors.
%     Samples = 5000 + 50;
%     error = zeros(Samples,1);
%     aux = weights(:,s+1);
%     weights = zeros(order, Samples);
%     weights(:,order) = aux;

%     signal_d_hat_50 = zeros(size(signal_d));
%     % LMS algorithm
%     for s = order:Samples
%         aux = signal_x(s:-1:s-order+1);
%         % Filtering the signal
%         signal_d_hat_50(s-order+1) = weights(:,s)'*aux;
%         % The equalizer does not know the original signal
%         error(s) = qammod(qamdemod(signal_x(s-order+1),QAM),QAM) - weights(:,s)'*aux;
%         % Recursive expression.
%         weights(:,s+1) = weights(:,s) + 2 * mi * conj(error(s)) * aux;
%     end

%     % Training with 150 Samples

%     % Training Phase

%     % Number of samples
%     Samples = 150;

%     % Defining the mse error and filter coeficients vectors.
%     error = zeros(Samples,1);
%     weights = zeros(order, Samples);

%     % Defining the energy of the noise vector.
%     QAM_train = 4;
%     signal_d_train = randi([0,QAM_train - 1],[Samples 1]);
%     signal_d_train = (1/sqrt(2)) * qammod(signal_d_train,QAM_train); 

%     % Convolving the channel and the signal.
%     Hz = [0.5 1.2 1.5 -1];
%     signal_x_train = filter(Hz,1,signal_d_train);

%     % Training noise
%     snr = 10^(inf/10);
%     energy_symbol = mean(abs(signal_x_train(:)).^2); % Energy symbol pilot. 
%     var_noise = energy_symbol .*  1/snr; % Variance of the noise.
%     noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

%     % Generating the noisy received signal.
%     signal_x_train = signal_x_train + noise;
%     % LMS algorithm
%     for s = order:Samples
%         aux = signal_x_train(s:-1:s-order+1);
%         error(s) = signal_d_train(s-order+1) - weights(:,s)'*aux;
%         % Recursive expression.
%         weights(:,s+1) = weights(:,s) + 2 * mi * conj(error(s)) * aux;
%     end

%     % Transmission Phase

%     % Defining the mse error and filter coeficients vectors.
%     Samples = 5000 + 50;
%     error = zeros(Samples,1);
%     aux = weights(:,s+1);
%     weights = zeros(order, Samples);
%     weights(:,order) = aux;

%     signal_d_hat_150 = zeros(size(signal_d));
%     % LMS algorithm
%     for s = order:Samples
%         aux = signal_x(s:-1:s-order+1);
%         % Filtering the signal
%         signal_d_hat_150(s-order+1) = weights(:,s)'*aux;
%         % The equalizer does not know the original signal
%         error(s) = qammod(qamdemod(signal_x(s-order+1),QAM),QAM) - weights(:,s)'*aux;
%         % Recursive expression.
%         weights(:,s+1) = weights(:,s) + 2 * mi * conj(error(s)) * aux;
%     end

%     % Training with 300 Samples

%     % Training Phase

%     % Number of samples
%     Samples = 300;

%     % Defining the mse error and filter coeficients vectors.
%     error = zeros(Samples,1);
%     weights = zeros(order, Samples);

%     % Defining the energy of the noise vector.
%     QAM_train = 4;
%     signal_d_train = randi([0,QAM_train - 1],[Samples 1]);
%     signal_d_train = (1/sqrt(2)) * qammod(signal_d_train,QAM_train); 

%     % Convolving the channel and the signal.
%     Hz = [0.5 1.2 1.5 -1];
%     signal_x_train = filter(Hz,1,signal_d_train);

%     % Training noise
%     snr = 10^(inf/10);
%     energy_symbol = mean(abs(signal_x_train(:)).^2); % Energy symbol pilot. 
%     var_noise = energy_symbol .*  1/snr; % Variance of the noise.
%     noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

%     % Generating the noisy received signal.
%     signal_x_train = signal_x_train + noise;
%     % LMS algorithm
%     for s = order:Samples
%         aux = signal_x_train(s:-1:s-order+1);
%         error(s) = signal_d_train(s-order+1) - weights(:,s)'*aux;
%         % Recursive expression.
%         weights(:,s+1) = weights(:,s) + 2 * mi * conj(error(s)) * aux;
%     end

%     % Transmission Phase

%     % Defining the mse error and filter coeficients vectors.
%     Samples = 5000 + 50;
%     error = zeros(Samples,1);
%     aux = weights(:,s+1);
%     weights = zeros(order, Samples);
%     weights(:,order) = aux;

%     signal_d_hat_300 = zeros(size(signal_d));
%     % LMS algorithm
%     for s = order:Samples
%         aux = signal_x(s:-1:s-order+1);
%         % Filtering the signal
%         signal_d_hat_300(s-order+1) = weights(:,s)'*aux;
%         % The equalizer does not know the original signal
%         error(s) = qammod(qamdemod(signal_x(s-order+1),QAM),QAM) - weights(:,s)'*aux;
%         % Recursive expression.
%         weights(:,s+1) = weights(:,s) + 2 * mi * conj(error(s)) * aux;
%     end

%     % Training with 500 Samples

%     % Training Phase

%     % Number of samples
%     Samples = 500;

%     % Defining the mse error and filter coeficients vectors.
%     error = zeros(Samples,1);
%     weights = zeros(order, Samples);

%     % Defining the energy of the noise vector.
%     QAM_train = 4;
%     signal_d_train = randi([0,QAM_train - 1],[Samples 1]);
%     signal_d_train = qammod(signal_d_train,QAM_train); 

%     % Convolving the channel and the signal.
%     Hz = [0.5 1.2 1.5 -1];
%     signal_x_train = filter(Hz,1,signal_d_train);

%     % Training noise
%     snr = 10^(inf/10);
%     energy_symbol = mean(abs(signal_x_train(:)).^2); % Energy symbol pilot. 
%     var_noise = energy_symbol .*  1/snr; % Variance of the noise.
%     noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

%     % Generating the noisy received signal.
%     signal_x_train = signal_x_train + noise;
%     % LMS algorithm
%     for s = order:Samples
%         aux = signal_x_train(s:-1:s-order+1);
%         error(s) = signal_d_train(s-order+1) - weights(:,s)'*aux;
%         % Recursive expression.
%         weights(:,s+1) = weights(:,s) + 2 * mi * conj(error(s)) * aux;
%     end

%     % Transmission Phase

%     % Defining the mse error and filter coeficients vectors.
%     Samples = 5000 + 50;
%     error = zeros(Samples,1);
%     aux = weights(:,s+1);
%     weights = zeros(order, Samples);
%     weights(:,order) = aux;

%     signal_d_hat_500 = zeros(size(signal_d));
%     % LMS algorithm
%     for s = order:Samples
%         aux = signal_x(s:-1:s-order+1);
%         % Filtering the signal
%         signal_d_hat_500(s-order+1) = weights(:,s)'*aux;
%         % The equalizer does not know the original signal
%         error(s) = qammod(qamdemod(signal_x(s-order+1),QAM),QAM) - weights(:,s)'*aux;
%         % Recursive expression.
%         weights(:,s+1) = weights(:,s) + 2 * mi * conj(error(s)) * aux;
%     end

%     % Temporal Evolution 
%     range_of_sampling = 4900:5000;
%     [~,~,delay] = alignsignals(qamdemod(signal_d,QAM),qamdemod(signal_d_hat_500,QAM));

%     aux = qamdemod(signal_d_hat_50,QAM);
%     aux = circshift(aux,delay);
%     aux_50 = aux(range_of_sampling);

%     aux = qamdemod(signal_d_hat_150,QAM);
%     aux = circshift(aux,delay);
%     aux_150 = aux(range_of_sampling);

%     aux = qamdemod(signal_d_hat_300,QAM);
%     aux = circshift(aux,delay);
%     aux_300 = aux(range_of_sampling);

%     aux = qamdemod(signal_d_hat_500,QAM);
%     aux = circshift(aux,delay);
%     aux_500 = aux(range_of_sampling);

%     figure
%     subplot(221)
%     txt = ['Original Signal'];
%     stem(range_of_sampling, qamdemod(signal_d(range_of_sampling),QAM),'-','color', [0.3010 0.7450 0.9330], "linewidth", 1, "markersize", 3, "DisplayName", txt);
%     hold on;
%     txt = ['Estimated Signal'];
%     stem(range_of_sampling, aux_50,'-','color', [0.4660 0.6740 0.1880], "linewidth", 1, "markersize", 3, "DisplayName", txt);
%     hold off;
%     title('50 Training Samples');
%     xlabel('Sample');
%     xlim([min(range_of_sampling) max(range_of_sampling)]);
%     ylabel('Magnitude');
%     legend_copy = legend("location", "southwest");
%     set (legend_copy, "fontsize", 6);
%     grid on;
%     subplot(222)
%     txt = ['Original Signal'];
%     stem(range_of_sampling, qamdemod(signal_d(range_of_sampling),QAM),'-','color', [0.3010 0.7450 0.9330], "linewidth", 1, "markersize", 3, "DisplayName", txt);
%     hold on;
%     txt = ['Estimated Signal'];
%     stem(range_of_sampling, aux_150,'-','color', [0.4660 0.6740 0.1880], "linewidth", 1, "markersize", 3, "DisplayName", txt);
%     hold off;
%     title('150 Training Samples');
%     xlabel('Sample');
%     xlim([min(range_of_sampling) max(range_of_sampling)]);
%     ylabel('Magnitude');
%     legend_copy = legend("location", "southwest");
%     set (legend_copy, "fontsize", 6);
%     grid on;
%     subplot(223)
%     txt = ['Original Signal'];
%     stem(range_of_sampling, qamdemod(signal_d(range_of_sampling),QAM),'-','color', [0.3010 0.7450 0.9330], "linewidth", 1, "markersize", 3, "DisplayName", txt);
%     hold on;
%     txt = ['Estimated Signal'];
%     stem(range_of_sampling, aux_300,'-','color', [0.4660 0.6740 0.1880], "linewidth", 1, "markersize", 3, "DisplayName", txt);
%     hold off;
%     title('300 Training Samples');
%     xlabel('Sample');
%     xlim([min(range_of_sampling) max(range_of_sampling)]);
%     ylabel('Magnitude');
%     legend_copy = legend("location", "southwest");
%     set (legend_copy, "fontsize", 6);
%     grid on;
%     subplot(224)
%     txt = ['Original Signal'];
%     stem(range_of_sampling, qamdemod(signal_d(range_of_sampling),QAM),'-','color', [0.3010 0.7450 0.9330], "linewidth", 1, "markersize", 3, "DisplayName", txt);
%     hold on;
%     txt = ['Estimated Signal'];
%     stem(range_of_sampling, aux_500,'-','color', [0.4660 0.6740 0.1880], "linewidth", 1, "markersize", 3, "DisplayName", txt);
%     hold off;
%     title('500 Training Samples');
%     xlabel('Sample');
%     xlim([min(range_of_sampling) max(range_of_sampling)]);
%     ylabel('Magnitude');
%     legend_copy = legend("location", "southwest");
%     set (legend_copy, "fontsize", 6);
%     grid on;
%     %saveas(gcf,'L3Q6_B_t.png')

%     % (c) --------------------------------------
%     disp('c')
%     % Training Stage  

%     % Learning rate
%     mi = 0.4;
%     gamma = 1e-3;
%     % Filter order
%     order = 15 + 1;
%     % Number of samples
%     Samples = 500;
%     % Defining the mse error and filter coeficients vectors.
%     error = zeros(Samples,1);
%     weights = zeros(order, Samples);

%     % Defining the energy of the noise vector.
%     SNR = 30;
%     QAM_train = 4;
%     signal_d_train = randi([0,QAM_train - 1],[Samples 1]); % The same pilot for every pilot frame and block.
%     signal_d_train = qammod(signal_d_train,QAM_train); % 4-QAM Pilot Signal.

%     % Convolving the channel and the signal.
%     Hz = [0.5 1.2 1.5 -1];
%     signal_x_train = filtfilt(Hz,1,signal_d_train);

%     % Training noise
%     snr = 10^(inf/10);
%     energy_symbol = mean(abs(signal_x_train(:)).^2); % Energy symbol pilot. 
%     var_noise = energy_symbol .*  1/snr; % Variance of the noise.
%     noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

%     % Generating the noisy received signal.
%     signal_x_train = signal_x_train + noise;
%     for s = order:Samples
%         aux = signal_x_train(s:-1:s-order+1);
%         mi_normalized = mi/(gamma + norm(aux)^2);
%         error(s) = signal_d_train(s-order+1) - weights(:,s)'*aux;
%         % Recursive expression.
%         weights(:,s+1) = weights(:,s) + mi_normalized * conj(error(s)) * aux;
%     end

%     % Transmission Stage 

%     % Number of samples
%     Samples = 5000 + 50;
%     % Defining the mse error and filter coeficients vectors.
%     error = zeros(Samples,1);
%     weights = zeros(order, Samples);

%     % Defining the energy of the noise vector.
%     SNR = 30;
%     QAM = 256;
%     signal_d = randi([0,QAM - 1],[Samples 1]); % The same pilot for every pilot frame and block.
%     signal_d = qammod(signal_d,QAM); % 4-QAM Pilot Signal.

%     % Convolving the channel and the signal.
%     Hz = [0.5 1.2 1.5 -1];
%     signal_x = filtfilt(Hz,1,signal_d);

%     % Training noise
%     snr = 10^(SNR/10);
%     energy_symbol = mean(abs(signal_x(:)).^2); % Energy symbol pilot. 
%     var_noise = energy_symbol .*  1/snr; % Variance of the noise.
%     noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

%     % Generating the noisy received signal.
%     signal_x = signal_x + noise;
%     signal_d_hat = zeros(size(signal_d));
%     % NLMS algorithm
%     for s = order:Samples
%         aux = signal_x(s:-1:s-order+1);
%         mi_normalized = mi/(gamma + norm(aux)^2);
%         % Filtering the signal
%         signal_d_hat(s-order+1) = weights(:,s)'*aux;
%         % The equalizer does know the original signal
%         %error(s) = signal_d(s-order+1) - weights(:,s)'*aux;
%         % The equalizer does not know the original signal
%         error(s) = qammod(qamdemod(signal_x(s-order+1),QAM),QAM) - weights(:,s)'*aux;
%         % Recursive expression.
%         weights(:,s+1) = weights(:,s) + mi_normalized * conj(error(s)) * aux;
%     end

%     % MSE Curve
%     figure
%     semilogy(1:Samples, abs(error).^2,'-','color', [0.3010 0.7450 0.9330], "linewidth", 1, "markersize", 8);
%     title('NLMS Behavior');
%     xlabel('Samples');
%     xlim([0 5000])
%     ylabel('MSE');
%     grid on;
%     % saveas(gcf,'L3Q6_C_mse.png')

%     % Temporal Evolution
%     aux = qamdemod(signal_d_hat,QAM);
%     aux1 = aux(1:100);
%     aux2 = aux(4900:5000);

%     figure
%     subplot(211)
%     txt = ['Original Signal'];
%     stem(1:100, qamdemod(signal_d(1:100),QAM),'-','color', [0.3010 0.7450 0.9330], "linewidth", 1, "markersize", 3, "DisplayName", txt);
%     hold on;
%     txt = ['Estimated Signal'];
%     stem(1:100, aux1,'-','color', [0.4660 0.6740 0.1880], "linewidth", 1, "markersize", 3, "DisplayName", txt);
%     hold off;
%     title('Temporal Evolution of First Samples');
%     xlabel('Sample');
%     xlim([0 100])
%     ylabel('Magnitude');
%     legend_copy = legend("location", "southwest");
%     set (legend_copy, "fontsize", 6);
%     grid on;
%     subplot(212)
%     txt = ['Original Signal'];
%     stem(4900:5000, qamdemod(signal_d(4900:5000),QAM),'-','color', [0.3010 0.7450 0.9330], "linewidth", 1, "markersize", 3, "DisplayName", txt);
%     hold on;
%     txt = ['Estimated Signal'];
%     stem(4900:5000, aux2,'-','color', [0.4660 0.6740 0.1880], "linewidth", 1, "markersize", 3, "DisplayName", txt);
%     hold off;
%     title('Temporal Evolution of Last Samples');
%     xlabel('Sample');
%     xlim([4900 5000])
%     ylabel('Magnitude');
%     legend_copy = legend("location", "southwest");
%     set (legend_copy, "fontsize", 6);
%     grid on;
%     % saveas(gcf,'L3Q6_C_t.png')

%     % Constellation
%     % https://www.mathworks.com/help/comm/gs/examine-16-qam-using-matlab.html
%     %figure
%     %subplot(221)
%     %plot(signal_d_train,'.','color', [0.3010 0.7450 0.9330],"markersize", 3)
%     %title('Training Signal');
%     %xlabel('In Phase');
%     %ylabel('Quadrature');
%     %grid on;
%     %subplot(222)
%     %plot(signal_d,'.','color', [0.3010 0.7450 0.9330],"markersize", 3)
%     %title('Original Signal');
%     %xlabel('In Phase');
%     %ylabel('Quadrature');
%     %grid on;
%     %subplot(223)
%     %plot(signal_x,'.','color', [0.3010 0.7450 0.9330],"markersize", 3)
%     %title('Transmitted Signal');
%     %xlabel('In Phase');
%     %ylabel('Quadrature');
%     %grid on;
%     %subplot(224)
%     %plot(signal_d_hat,'.','color', [0.3010 0.7450 0.9330],"markersize", 3)
%     %title('Filtered Signal');
%     %xlabel('In Phase');
%     %ylabel('Quadrature');
%     %grid on;
%     %saveas(gcf,'L3Q6_C_t.png')

%     % (d) --------------------------------------
%     disp('d')
%     % Simulation parameters
%     mi = 0.4;
%     gamma = 1e-3;
%     runs = 1; % runs = 1000;
%     QAM_train = 4;
%     snrs = [0 10 20 30];
%     % Filter order
%     order = 15 + 1;

%     tic;
%     % 4QAM
%     QAM = 4
%     SER = zeros(length(snrs),1);
%     for ii = 1:length(snrs)
%         snrs(ii)
%         for rr = 1:runs
%             % Training Stage  
%             % Number of samples
%             Samples = 500;
%             % Defining the mse error and filter coeficients vectors.
%             error = zeros(Samples,1);
%             weights = zeros(order, Samples);

%             % Defining the energy of the noise vector.
%             signal_d_train = randi([0,QAM_train - 1],[Samples 1]);
%             signal_d_train = qammod(signal_d_train,QAM_train);

%             % Convolving the channel and the signal.
%             Hz = [0.5 1.2 1.5 -1];
%             signal_x_train = filtfilt(Hz,1,signal_d_train);

%             % Training noise
%             snr = 10^(inf/10);
%             energy_symbol = mean(abs(signal_x_train(:)).^2); % Energy symbol pilot. 
%             var_noise = energy_symbol .*  1/snr; % Variance of the noise.
%             noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

%             % Generating the noisy received signal.
%             signal_x_train = signal_x_train + noise;

%             for s = order:Samples
%                 aux = signal_x_train(s:-1:s-order+1);
%                 mi_normalized = mi/(gamma + norm(aux)^2);
%                 error(s) = signal_d_train(s-order+1) - weights(:,s)'*aux;
%                 % Recursive expression.
%                 weights(:,s+1) = weights(:,s) + mi_normalized * conj(error(s)) * aux;
%             end

%             % Transmission Stage 
%             % Number of samples
%             Samples = 5000;
%             % Defining the mse error and filter coeficients vectors.
%             error = zeros(Samples,1);
%             weights = zeros(order, Samples);

%             % Defining the energy of the noise vector.
%             signal_d = randi([0,QAM - 1],[Samples 1]);
%             signal_d = qammod(signal_d,QAM);

%             % Convolving the channel and the signal.
%             Hz = [0.5 1.2 1.5 -1];
%             signal_x = filtfilt(Hz,1,signal_d);

%             % Transmission noise
%             snr = 10^(snrs(ii)/10);
%             energy_symbol = mean(abs(signal_x(:)).^2); % Energy symbol pilot. 
%             var_noise = energy_symbol .*  1/snr; % Variance of the noise.
%             noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

%             % Generating the noisy received signal.
%             signal_x = signal_x + noise;
%             signal_d_hat = zeros(size(signal_d));
%             % NLMS algorithm
%             for s = order:Samples
%                 aux = signal_x(s:-1:s-order+1);
%                 mi_normalized = mi/(gamma + norm(aux)^2);
%                 % Filtering the signal
%                 signal_d_hat(s-order+1) = weights(:,s)'*aux;
%                 % The equalizer does not know the original signal
%                 error(s) = qammod(qamdemod(signal_x(s-order+1),QAM),QAM) - weights(:,s)'*aux;
%                 % Recursive expression.
%                 weights(:,s+1) = weights(:,s) + mi_normalized * conj(error(s)) * aux;
%             end
%             aux1 = qamdemod(signal_d,QAM);
%             aux2 = qamdemod(signal_d_hat,QAM);
%             SER(ii,1) = SER(ii,1) + sum(aux1~=aux2)/length(aux1);
%         end
%     end
%     SER = SER/runs;
%     figure
%     txt = ['4QAM Signal'];
%     semilogy(snrs, SER,'-','color', [0.3010 0.7450 0.9330], "linewidth", 3, "markersize", 8, "DisplayName", txt);
%     hold on;

%     % 16QAM 
%     QAM = 16
%     SER = zeros(length(snrs),1);
%     for ii = 1:length(snrs)
%         snrs(ii)
%         for rr = 1:runs
%             % Training Stage  
%             % Number of samples
%             Samples = 500;
%             % Defining the mse error and filter coeficients vectors.
%             error = zeros(Samples,1);
%             weights = zeros(order, Samples);

%             % Defining the energy of the noise vector.
%             signal_d_train = randi([0,QAM_train - 1],[Samples 1]);
%             signal_d_train = qammod(signal_d_train,QAM_train);

%             % Convolving the channel and the signal.
%             Hz = [0.5 1.2 1.5 -1];
%             signal_x_train = filtfilt(Hz,1,signal_d_train);

%             % Training noise
%             snr = 10^(inf/10);
%             energy_symbol = mean(abs(signal_x_train(:)).^2); % Energy symbol pilot. 
%             var_noise = energy_symbol .*  1/snr; % Variance of the noise.
%             noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

%             % Generating the noisy received signal.
%             signal_x_train = signal_x_train + noise;

%             for s = order:Samples
%                 aux = signal_x_train(s:-1:s-order+1);
%                 mi_normalized = mi/(gamma + norm(aux)^2);
%                 error(s) = signal_d_train(s-order+1) - weights(:,s)'*aux;
%                 % Recursive expression.
%                 weights(:,s+1) = weights(:,s) + mi_normalized * conj(error(s)) * aux;
%             end

%             % Transmission Stage 
%             % Number of samples
%             Samples = 5000;
%             % Defining the mse error and filter coeficients vectors.
%             error = zeros(Samples,1);
%             weights = zeros(order, Samples);

%             % Defining the energy of the noise vector.
%             signal_d = randi([0,QAM - 1],[Samples 1]);
%             signal_d = qammod(signal_d,QAM);

%             % Convolving the channel and the signal.
%             Hz = [0.5 1.2 1.5 -1];
%             signal_x = filtfilt(Hz,1,signal_d);

%             % Transmission noise
%             snr = 10^(snrs(ii)/10);
%             energy_symbol = mean(abs(signal_x(:)).^2); % Energy symbol pilot. 
%             var_noise = energy_symbol .*  1/snr; % Variance of the noise.
%             noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

%             % Generating the noisy received signal.
%             signal_x = signal_x + noise;
%             signal_d_hat = zeros(size(signal_d));
%             % NLMS algorithm
%             for s = order:Samples
%                 aux = signal_x(s:-1:s-order+1);
%                 mi_normalized = mi/(gamma + norm(aux)^2);
%                 % Filtering the signal
%                 signal_d_hat(s-order+1) = weights(:,s)'*aux;
%                 % The equalizer does not know the original signal
%                 error(s) = qammod(qamdemod(signal_x(s-order+1),QAM),QAM) - weights(:,s)'*aux;
%                 % Recursive expression.
%                 weights(:,s+1) = weights(:,s) + mi_normalized * conj(error(s)) * aux;
%             end
%             aux1 = qamdemod(signal_d,QAM);
%             aux2 = qamdemod(signal_d_hat,QAM);
%             SER(ii,1) = SER(ii,1) + sum(aux1~=aux2)/length(aux1);
%         end
%     end
%     SER = SER/runs;
%     txt = ['16QAM Signal'];
%     semilogy(snrs, SER,'-','color', [0 0.4470 0.7410], "linewidth", 3, "markersize", 8, "DisplayName", txt);
%     hold on;

%     % 64QAM 
%     QAM = 64
%     SER = zeros(length(snrs),1);
%     for ii = 1:length(snrs)
%         snrs(ii)
%         for rr = 1:runs
%             % Training Stage  
%             % Number of samples
%             Samples = 500;
%             % Defining the mse error and filter coeficients vectors.
%             error = zeros(Samples,1);
%             weights = zeros(order, Samples);

%             % Defining the energy of the noise vector.
%             signal_d_train = randi([0,QAM_train - 1],[Samples 1]);
%             signal_d_train = qammod(signal_d_train,QAM_train);

%             % Convolving the channel and the signal.
%             Hz = [0.5 1.2 1.5 -1];
%             signal_x_train = filtfilt(Hz,1,signal_d_train);

%             % Training noise
%             snr = 10^(inf/10);
%             energy_symbol = mean(abs(signal_x_train(:)).^2); % Energy symbol pilot. 
%             var_noise = energy_symbol .*  1/snr; % Variance of the noise.
%             noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

%             % Generating the noisy received signal.
%             signal_x_train = signal_x_train + noise;

%             for s = order:Samples
%                 aux = signal_x_train(s:-1:s-order+1);
%                 mi_normalized = mi/(gamma + norm(aux)^2);
%                 error(s) = signal_d_train(s-order+1) - weights(:,s)'*aux;
%                 % Recursive expression.
%                 weights(:,s+1) = weights(:,s) + mi_normalized * conj(error(s)) * aux;
%             end

%             % Transmission Stage 
%             % Number of samples
%             Samples = 5000;
%             % Defining the mse error and filter coeficients vectors.
%             error = zeros(Samples,1);
%             weights = zeros(order, Samples);

%             % Defining the energy of the noise vector.
%             signal_d = randi([0,QAM - 1],[Samples 1]);
%             signal_d = qammod(signal_d,QAM);

%             % Convolving the channel and the signal.
%             Hz = [0.5 1.2 1.5 -1];
%             signal_x = filtfilt(Hz,1,signal_d);

%             % Transmission noise
%             snr = 10^(snrs(ii)/10);
%             energy_symbol = mean(abs(signal_x(:)).^2); % Energy symbol pilot. 
%             var_noise = energy_symbol .*  1/snr; % Variance of the noise.
%             noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));
%             % Generating the noisy received signal.
%             signal_x = signal_x + noise;
%             signal_d_hat = zeros(size(signal_d));
%             % NLMS algorithm
%             for s = order:Samples
%                 aux = signal_x(s:-1:s-order+1);
%                 mi_normalized = mi/(gamma + norm(aux)^2);
%                 % Filtering the signal
%                 signal_d_hat(s-order+1) = weights(:,s)'*aux;
%                 % The equalizer does not know the original signal
%                 error(s) = qammod(qamdemod(signal_x(s-order+1),QAM),QAM) - weights(:,s)'*aux;
%                 % Recursive expression.
%                 weights(:,s+1) = weights(:,s) + mi_normalized * conj(error(s)) * aux;
%             end
%             aux1 = qamdemod(signal_d,QAM);
%             aux2 = qamdemod(signal_d_hat,QAM);
%             SER(ii,1) = SER(ii,1) + sum(aux1~=aux2)/length(aux1);
%         end
%     end
%     SER = SER/runs;
%     txt = ['64QAM Signal'];
%     semilogy(snrs, SER,'-','color', [0.8500 0.3250 0.0980], "linewidth", 3, "markersize", 8, "DisplayName", txt);
%     hold on;

%     % 256QAM 
%     QAM = 256
%     SER = zeros(length(snrs),1);
%     for ii = 1:length(snrs)
%         snrs(ii)
%         for rr = 1:runs
%             % Training Stage  
%             % Number of samples
%             Samples = 500;
%             % Defining the mse error and filter coeficients vectors.
%             error = zeros(Samples,1);
%             weights = zeros(order, Samples);

%             % Defining the energy of the noise vector.
%             signal_d_train = randi([0,QAM_train - 1],[Samples 1]);
%             signal_d_train = qammod(signal_d_train,QAM_train);

%             % Convolving the channel and the signal.
%             Hz = [0.5 1.2 1.5 -1];
%             signal_x_train = filtfilt(Hz,1,signal_d_train);

%             % Training noise
%             snr = 10^(inf/10);
%             energy_symbol = mean(abs(signal_x_train(:)).^2); % Energy symbol pilot. 
%             var_noise = energy_symbol .*  1/snr; % Variance of the noise.
%             noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

%             % Generating the noisy received signal.
%             signal_x_train = signal_x_train + noise;

%             for s = order:Samples
%                 aux = signal_x_train(s:-1:s-order+1);
%                 mi_normalized = mi/(gamma + norm(aux)^2);
%                 error(s) = signal_d_train(s-order+1) - weights(:,s)'*aux;
%                 % Recursive expression.
%                 weights(:,s+1) = weights(:,s) + mi_normalized * conj(error(s)) * aux;
%             end

%             % Transmission Stage 
%             % Number of samples
%             Samples = 5000;
%             % Defining the mse error and filter coeficients vectors.
%             error = zeros(Samples,1);
%             weights = zeros(order, Samples);

%             % Defining the energy of the noise vector.
%             signal_d = randi([0,QAM - 1],[Samples 1]);
%             signal_d = qammod(signal_d,QAM);

%             % Convolving the channel and the signal.
%             Hz = [0.5 1.2 1.5 -1];
%             signal_x = filtfilt(Hz,1,signal_d);

%             % Transmission noise
%             snr = 10^(snrs(ii)/10);
%             energy_symbol = mean(abs(signal_x(:)).^2); % Energy symbol pilot. 
%             var_noise = energy_symbol .*  1/snr; % Variance of the noise.
%             noise = sqrt(var_noise/2) * (randn(Samples,1) + 1i*randn(Samples,1));

%             % Generating the noisy received signal.
%             signal_x = signal_x + noise;
%             signal_d_hat = zeros(size(signal_d));
%             % NLMS algorithm
%             for s = order:Samples
%                 aux = signal_x(s:-1:s-order+1);
%                 mi_normalized = mi/(gamma + norm(aux)^2);
%                 % Filtering the signal
%                 signal_d_hat(s-order+1) = weights(:,s)'*aux;
%                 % The equalizer does not know the original signal
%                 error(s) = qammod(qamdemod(signal_x(s-order+1),QAM),QAM) - weights(:,s)'*aux;
%                 % Recursive expression.
%                 weights(:,s+1) = weights(:,s) + mi_normalized * conj(error(s)) * aux;
%             end
%             aux1 = qamdemod(signal_d,QAM);
%             aux2 = qamdemod(signal_d_hat,QAM);
%             SER(ii,1) = SER(ii,1) + sum(aux1~=aux2)/length(aux1);
%         end
%     end

%     t = toc; 
    
%     disp(t)

%     SER = SER/runs;
%     txt = ['256QAM Signal'];
%     semilogy(snrs, SER,'-','color', [0.4660 0.6740 0.1880], "linewidth", 3, "markersize", 8, "DisplayName", txt);
%     hold off;
%     legend_copy = legend("location", "southwest");
%     set (legend_copy, "fontsize", 12);
%     grid on;
%     title('SER vs. SNR for different constelattions');
%     xlabel('SNR (dB)');
%     ylabel('SER');
%     % saveas(gcf,'L3Q6_D_ser.png')
    
% end


% %% HOMEWORK 4 - PROBLEM 1
% function hw4p1(varargin)
    
%     % Forgeting rate
%     lambda = 0.98;
%     delta = 1;
%     % Filter order
%     order = 2 + 1;


%     % Number of samples
%     Samples = 100;
%     % Defining the mse error and filter coeficients vectors.
%     error = zeros(Samples,1);
%     weights = zeros(order, Samples);
%     weights(1,1) = 1; 
%     % Generating the sinoide signal.
%     t = linspace(-pi,pi,Samples).';
%     signal_d = cos(pi*t/3);
%     % Defining the energy of the noise vector.
%     SNR_dB = 15;
%     SNR_li = 10^(SNR_dB/10);
%     variance_noise = 1/SNR_li;
%     noise = sqrt(variance_noise/2).*randn(Samples,1);
%     % Generating the noisy received signal.
%     signal_x = signal_d + noise;

%     % Deterministic correlation matrix initialization
%     y = zeros(Samples,1);
%     Rd = delta*eye(order); 
%     for ss = 2:(Samples - order - 1)
%         % Deterministic correlation matrix inverse
%         Rd = (1/lambda)*(Rd - (Rd*signal_x(ss:ss+order-1)*signal_x(ss:ss+order-1)'*Rd)/(lambda + signal_x(ss:ss+order-1)'*Rd*signal_x(ss:ss+order-1)));
%         % Error between the desired signal and the filtered signal.
%         error(ss) = signal_d(ss) - weights(:,ss-1)' * signal_x(ss:ss+order-1); 
%         % Recursive expression.
%         weights(:,ss) = weights(:,ss-1) + Rd*error(ss)*signal_x(ss:ss+order-1);
%         weights(1,ss) = 1;
%         % Output
%         y(ss) = weights(:,ss-1)' * signal_x(ss:ss+order-1);
%     end

%     aux = weights(:,1:10).';
%     tabela = [aux(1,:);aux(2,:);aux(3,:);aux(4,:);aux(5,:);aux(6,:);aux(7,:);aux(8,:);aux(9,:);aux(10,:);];
%     Tabela = array2table(tabela);
%     Tabela.Properties.RowNames = {'i = 1' 'i = 2' 'i = 3' 'i = 4' 'i = 5' 'i = 6' 'i = 7' 'i = 8' 'i = 9' 'i = 10'};
%     Tabela.Properties.VariableNames = {'1st Coeff' '2nd Coeff' '3rd Coeff'};
%     % table2latex(Tabela,'tabela.tex');

%     % MSE Curve
%     figure
%     txt = ['Filtered Signal'];
%     plot(y,'-','color', [0.4660 0.6740 0.1880], "linewidth", 2, "markersize", 8, "DisplayName", txt);
%     hold on;
%     txt = ['Source Signal'];
%     plot(signal_d,'-','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName", txt);
%     hold off;
%     title('Predicted Signal vs. Source Signal');
%     xlabel('Samples');
%     ylabel('Magnitude');
%     legend_copy = legend("location", "northeast");
%     set (legend_copy, "fontsize", 12);
%     grid on;
%     % saveas(gcf,'L4Q1.png')

% end


% %% HOMEWORK 4 - PROBLEM 3
% function hw4p3(varargin)
%     % Forgeting rate
%     lambda = 0.99;
%     % Filter order
%     order = 3;
%     % SNR (dB) 
%     SNR_dB = 3;

%     % Number of samples
%     Samples = 1010;
%     % Defining the mse error and filter coeficients vectors.
%     error = zeros(Samples,1);
%     weights = zeros(order, Samples);
%     % Defining the energy of the noise vector.
%     SNR_li = 10^(SNR_dB/10);
%     variance_noise = 1/SNR_li;
%     noise = sqrt(variance_noise/2).*randn(Samples,1);
%     % Generating the sinoide signal.
%     t = linspace(-pi,pi,Samples).';
%     signal_d = sin(2*pi*t);
%     % Generating the noisy received signal.
%     signal_x = signal_d + noise;
%     % Defining delta by the inverse of the signal energy
%     delta  = 1/(sum(signal_x.^2)/length(signal_x));

%     % Deterministic correlation matrix initialization
%     Rd = delta*eye(order); 
%     %signal_d = signal_d(order:end,1); 
%     for ss = 2:(Samples - order - 1)
%         % Deterministic correlation matrix inverse
%         Rd = (1/lambda)*(Rd - (Rd*signal_x(ss:ss+order-1)*signal_x(ss:ss+order-1)'*Rd)/(lambda + signal_x(ss:ss+order-1)'*Rd*signal_x(ss:ss+order-1)));
%         % Error between the desired signal and the filtered signal.
%         error(ss) = signal_d(ss) - weights(:,ss-1)' * signal_x(ss:ss+order-1); 
%         % Recursive expression.
%         weights(:,ss) = weights(:,ss-1) + Rd*error(ss)*signal_x(ss:ss+order-1);
%     end
%     weights = flip(weights); 
%     % I first implemented thinking of python notation, later I found out that
%     % the reference book defines the order a bit different from what I usually
%     % work. So to make the code close to Diniz notation the 'order + 1' is
%     % needed.
%     % Coefficients Curve
%     figure
%     subplot(2,1,1)
%     txt = ['w0'];
%     plot(1:Samples, weights(1,:),'-','color', [0.3010 0.7450 0.9330], "linewidth", 2, "markersize", 8, "DisplayName",txt);
%     hold on;
%     txt = ['w1'];
%     plot(1:Samples, weights(2,:),'-','color', [0.4660 0.6740 0.1880], "linewidth", 2, "markersize", 8, "DisplayName",txt);
%     hold off;
%     title('RLS Coefficients Behavior');
%     xlabel('Samples');
%     xlim([0 1000]);
%     ylabel('Magnitude');
%     legend_copy = legend("location", "northeast");
%     set (legend_copy, "fontsize", 8);
%     grid on;

%     % MSE Curve
%     subplot(2,1,2)
%     semilogy(1:Samples, error.^2,'-','color', [0.3010 0.7450 0.9330], "linewidth", 1, "markersize", 8);
%     title('RLS Behavior');
%     xlabel('Samples');
%     xlim([0 1000]);
%     ylabel('MSE');
%     grid on;
    %order_snr_forgeting
    % saveas(gcf,'L4Q3_rls_mse_3_3_99.png')
    
end


%% VERBOSE DETAILS
function export_fig(Activate, h, filename)
    if Activate
        savefig_tight(h, filename, 'both');
        filter_hw.verbose_save(filename);
    else
        pause(2)
        close(h);
    end
end

function verbose_save(filename)
    fprintf('Saving Results for:\n\t %s \n', filename);
end


%% SAVE DATA TO TXT FILE
function mat2txt(filename, X, permission, header)
% ND.MAT2TXT  Write a matrix X into a txt file
%   mat2txt(filename, X, 'w', header) - Overwite the file
%   mat2txt(filename, X, 'a', header) - Append to the file end
%
%   See also.
        [I, J] = size(X);
        fileID = fopen(filename, permission);
        fprintf(fileID, [repelem('-', strlength(header)+3), '\n', header, ...
                '\n', repelem('-', strlength(header)+3), '\n']);
        fprintf(fileID, 'X(%d, %d)\n', I, J);
            for ii = 1:I
                for jj = 1:J
                    fprintf(fileID, ' %2.0f', X(ii,jj));
                end
                fprintf(fileID, ';\n');
            end
        fprintf(fileID, '\n');
        fclose(fileID);
end



% end methods list


end
end