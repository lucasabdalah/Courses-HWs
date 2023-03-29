%% [TIP7200 - Processamento Digital de Sinais] - Backward Difference
% Authour: Lucas Abdalah
% ----------------------------
% 
% ex2_backwardDifference.m
% 2023/03/29 - v1
% 

%% Prepare Ambient 
clearvars;
close all;
% pause on


%% Synthetic signal Input
x = [1, 2, -1, 2, 1, 1, -2, 1];
y = backwardDifference(x);

figure(1);
stem(x, 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5);
hold on;
plot(y, ...
  'Color', 'red',...        
  'LineStyle', '-',...
  'LineWidth', 1.5,...
  'Marker', 'v',...
  'MarkerFaceColor', 'red',...
  'MarkerSize', 5);
hold off;
legend("x", "y", 'Location', 'northeastoutside','Orientation', 'Horizontal','Position', [0.5 0.47 0.0 1], 'Units','normalized');
legend boxoff
xlabel("n");
ylabel("Amplitude");
grid on;


%% Audio Signal Input 1
[x_fala, Fs_fala] = audioread("fala_sino.wav");
N_fala = length(x_fala);
t_fala = linspace(0, N_fala/Fs_fala, N_fala);
y_fala = backwardDifference(x_fala);

figure(2);
plot(t_fala, x_fala, 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.0);
hold on;
plot(t_fala, y_fala, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1.0);
hold off;
legend("x", "y", 'Location', 'northeastoutside','Orientation', 'Horizontal','Position', [0.5 0.47 0.0 1], 'Units','normalized');
legend boxoff
xlabel("Time (s)")
ylabel("Amplitude")
grid on;
axis tight


%% Audio Signal Input 2
[x_cantina, Fs_cantina] = audioread("cantinaband.wav");
N_cantina = length(x_cantina);
t_cantina = linspace(0, N_cantina/Fs_cantina, N_cantina);
y_cantina = backwardDifference(x_cantina);

figure(3);
plot(t_cantina, x_cantina, 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.0);
hold on;
plot(t_cantina, y_cantina, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 1.0);
hold off;
legend("x", "y", 'Location', 'northeastoutside','Orientation', 'Horizontal','Position', [0.5 0.47 0.0 1], 'Units','normalized');
legend boxoff
xlabel("Time (s)")
ylabel("Amplitude")
grid on;
axis tight


%% Save results
savefig_tight(figure(1), '../figures/ex2_synthetic', 'both');
savefig_tight(figure(2), '../figures/ex2_fala_sino', 'both');
savefig_tight(figure(3), '../figures/ex2_cantina', 'both');


%% Local Functions
function y = backwardDifference(x)
% backwardDifference - Simple backward difference
% y[n] = x[n] - x[n-1]
%
% backwardDifference(x) This function computes the backward difference 
% of a time series x. We consider the first sample for input and  output
% equal to the, i.e, y[1] = x[1]  
% 
% INPUTS: 
%       x : signal
%
% OUTPUTS:
%       y : backward difference output
%
  y = zeros(1, length(x));

  y(1) = x(1); % Initialize the output sequence with the first sample

  for n = 2:length(x)
      y(n) = x(n) - x(n-1); % Compute the backward difference for each sample
  end
end


%% References
% [1] A. V. Oppenheim, R. W. Schafer, and J. R. Buck, Discrete-Time Signal Processing, 3rd. ed, Prentice Hall Press, 2009.