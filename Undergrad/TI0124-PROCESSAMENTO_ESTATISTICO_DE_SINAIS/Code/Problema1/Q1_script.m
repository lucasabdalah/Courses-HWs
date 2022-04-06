%% AP2 de PES
% Questao 1
%
% Q1_script.m
%
% 2021/08/24 - Lucas Abdalah

close all; clearvars; clc; % Clear the matlab ambient

% To reproduce the same results
rng('default');

%% General setup
  L = 30; % Filter Length (Number of coefficients) 
  N = 5e2 + L; % Number of samples
  lim_z = [0, N/2];
  n = linspace(lim_z(1), lim_z(2),N);

%% Input signal: z(n)
  sigma_v2 = 1e-4; % Noise variance 
  theta = pi; % Phase constant 
  z = sin((2 * pi / 60) * n + theta); % Input signal: z(n)

%% Input signal with additive white Gaussian noise (AWGN): d(n) = z(n) + v(n)
  v = randn(1, N)*sqrt(sigma_v2);
  d = z + v;
  d_window = d(L + 1:end);

%% lms_regular
  mu = 8e-3; % Algorithm Step (\mu)
  % start_state = 'random';
  start_state = 'random';
  [y,e,w] = lms_regular(d, L, mu, start_state);

%% Display Plots: 3 subplots
  h = figure();
  subplot(3,1,1);
    plot((1:N-L), w);
    str = ['Filter coefficients: $\mu$ = ', num2str(mu)]; title(str,'interpreter','latex');
    str = ['Weights, $w_{n}$']; ylabel(str,'interpreter','latex');
    grid on

  subplot(3,1,2);
    plot((1:N-L), e.^2,...
        'Color', 'b',...
        'LineWidth', 0.5,...
        'LineStyle', '-');
    hold on
    plot((1:N-L), mean(e.^2) .* ones(size(e)),...
        'Color', 'r',...
        'LineWidth', 2.0,...
        'LineStyle', '--');
    set(gca, 'YScale', 'log');
    str = ['Squared Error  $(e^2)$: Mean = ', num2str(mean(e.^2))]; title(str,'interpreter','latex');
    str = {'$e^2(n)$','Mean'}; legend(str,'interpreter','latex')
    str = ['$e^2(n)$']; ylabel(str,'interpreter','latex');
    grid on

  subplot(3,1,3)
    plot((1:N-L), y, ...
        'Color', 'r',...
        'LineWidth', 2.0,...
        'LineStyle', '-');
    hold on;
    plot((1:N-L), d_window, ...
        'Color', 'b', ...
        'LineWidth', 2.0,...
        'LineStyle', ':');
    str = ['$d(n)$ vs $y$(n)']; title(str,'interpreter','latex');
    str = {'$y$','d'}; legend(str,'interpreter','latex')
    str = ['Samples, $n$']; xlabel(str,'interpreter','latex');
    str = ['$z(n)$']; ylabel(str,'interpreter','latex');
    grid on

%% Save figure
  saveas(h,'random_start.svg');

pause(10)
close all