%% [TIP7200 - Processamento Digital de Sinais] - Moving Average
% Authour: Lucas Abdalah
% ----------------------------
% 
% ex1_moving_average.m
% 2023/03/24 - v1
% 

%% Prepare Ambient 
clearvars;
close all; 
pause off


%% Synthetic signal Input
x = [1, 2, -1, 2, 1, 1, -2, 1];
L = 5;
shifted = (1:length(x)) + floor(L/2);
y = ma(x, L);
y_matlab = movmean(x, L);

figure(1);
subplot(2,1,1); % Signal vs. implemented MA  vs. Matlab MA
stem(x, 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5);
hold on;
plot(y, ...
  'Color', 'red',...        
  'LineStyle', '-',...
  'LineWidth', 1.5,...
  'Marker', 'v',...
  'MarkerFaceColor', 'red',...
  'MarkerSize', 5);
plot(y_matlab,...
  'Color', 'blue',...        
  'LineStyle', '--',...
  'LineWidth', 1.5,...
  'Marker', 'square',...
  'MarkerFaceColor', 'blue',...
  'MarkerSize', 5);
hold off;
grid on;
xlabel("n")
ylabel("Amplitude")
xlim([-1, 11])
legend("x", "y", "y\_matlab", 'Location', 'northeastoutside','Orientation', 'Horizontal','Position', [0.5 0.47 0.0 1], 'Units','normalized');
legend boxoff

subplot(2,1,2); % Signal vs. shifted implemented MA  vs. Matlab MA
stem(x, 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5);
hold on;
plot(y, ...
  'Color', 'red',...        
  'LineStyle', '-',...
  'LineWidth', 1.5,...
  'Marker', 'v',...
  'MarkerFaceColor', 'red',...
  'MarkerSize', 5);
plot(shifted, y_matlab,...
  'Color', 'blue',...        
  'LineStyle', '--',...
  'LineWidth', 1.5,...
  'Marker', 'square',...
  'MarkerFaceColor', 'blue',...
  'MarkerSize', 5);
hold off;
grid on;
xlim([-1, 11]);
xlabel("n");
ylabel("Amplitude");
annotation('textarrow',[0.48 0.42], [0.18 0.32],'String',"Sample shift [+2]");
annotation('rectangle',[.25 .30 .1 .13],'Color','red');
annotation('rectangle',[.76 .24 .1 .13],'Color','blue');


%% Audio Signal Input 1
[x_fala, Fs_fala] = audioread("fala_sino.wav");
N_fala = length(x_fala);
t_fala = linspace(0, N_fala/Fs_fala, N_fala);
Lw_fala= [5, 10, 50];
y_fala = zeros(length(Lw_fala), N_fala);

for i = 1:length(Lw_fala)
  disp(i);
  y_fala(i,:) = ma(x_fala, Lw_fala(i));
  %sound(y_fala(i,:), Fs_fala);
  pause();
end

figure(2);
subplot(3,1,1); 
plot(t_fala, x_fala, 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.0);
hold on;
plot(t_fala, y_fala(1,:), 'Color', 'red', 'LineStyle', '-', 'LineWidth', 1.0);
text(4.7, 0.6, sprintf("$L = %g$", Lw_fala(1)), "HorizontalAlignment", "center", "interpreter", "Latex");
hold off;
grid on;
ylabel("Amplitude")
legend("x", "y", 'Location', 'northeastoutside','Orientation', 'Horizontal','Position', [0.5 0.47 0.0 1], 'Units','normalized');
legend boxoff

subplot(3,1,2); 
plot(t_fala, x_fala,...
  'Color', 'black',...        
  'LineStyle', ':',...
  'LineWidth', 1.0);
hold on;
plot(t_fala, y_fala(2,:), ...
  'Color', 'red',...        
  'LineStyle', '-',...
  'LineWidth', 1.0);
text(4.7, 0.6, sprintf("$L = %g$", Lw_fala(2)), "HorizontalAlignment", "center", "interpreter", "Latex");
ylabel("Amplitude")
hold off;
grid on;

subplot(3,1,3); 
plot(t_fala, x_fala,...
  'Color', 'black',...        
  'LineStyle', ':',...
  'LineWidth', 1.0);
hold on;
plot(t_fala, y_fala(3,:), ...
  'Color', 'red',...        
  'LineStyle', '-',...
  'LineWidth', 1.0);
text(4.7, 0.6, sprintf("$L = %g$", Lw_fala(3)), "HorizontalAlignment", "center", "interpreter", "Latex");
ylabel("Amplitude")
hold off;
xlabel("Time (s)")
grid on;


%% Audio Signal Input 2
[x_cantina, Fs_cantina] = audioread("cantinaband.wav");
N_cantina = length(x_cantina);
t_cantina = linspace(0, N_cantina/Fs_cantina, N_cantina);
Lw_cantina= [5, 10, 50];
y_cantina = zeros(length(Lw_cantina), N_cantina);

for i = 1:length(Lw_cantina)
  disp(i)
  y_cantina(i,:) = ma(x_cantina, Lw_cantina(i));
  %sound(y_cantina(i,:), Fs_cantina)
  pause()
end

figure(3);
subplot(3,1,1); 
plot(t_cantina, x_cantina,...
  'Color', 'black',...        
  'LineStyle', ':',...
  'LineWidth', 1.0);
hold on;
plot(t_cantina, y_cantina(1,:), ...
  'Color', 'red',...        
  'LineStyle', '-',...
  'LineWidth', 1.0)
hold off;
grid on;
text(0.25, 0.15, sprintf("$L = %g$", Lw_cantina(1)), "HorizontalAlignment", "center", "interpreter", "Latex");
ylabel("Amplitude")
legend("x", "y", 'Location', 'northeastoutside','Orientation', 'Horizontal','Position', [0.5 0.47 0.0 1], 'Units','normalized');
legend boxoff
subplot(3,1,2); 
plot(t_cantina, x_cantina,...
  'Color', 'black',...        
  'LineStyle', ':',...
  'LineWidth', 1.0);
hold on;
plot(t_cantina, y_cantina(2,:), ...
  'Color', 'red',...        
  'LineStyle', '-',...
  'LineWidth', 1.0);
hold off;
grid on;
ylabel("Amplitude")
text(0.25, 0.15, sprintf("$L = %g$", Lw_cantina(2)), "HorizontalAlignment", "center", "interpreter", "Latex");
subplot(3,1,3); 
plot(t_cantina, x_cantina,...
  'Color', 'black',...        
  'LineStyle', ':',...
  'LineWidth', 1.0);
hold on;
plot(t_cantina, y_cantina(3,:), ...
  'Color', 'red',...        
  'LineStyle', '-',...
  'LineWidth', 1.0);
hold off;
grid on;
xlabel("Time (s)")
ylabel("Amplitude")
text(0.25, 0.15, sprintf("$L = %g$", Lw_cantina(3)), "HorizontalAlignment", "center", "interpreter", "Latex");


%% Save results
savefig_tight(figure(1), '../figures/ex1_synthetic', 'both');
savefig_tight(figure(2), '../figures/ex1_fala_sino', 'both');
savefig_tight(figure(3), '../figures/ex1_cantina', 'both');


%% Local Functions
function y = ma(x, L)
% ma - Simple Moving Average
%
% ma(x, L) This function computes the moving average of a 
% time series x with a specified window size. Time series length must 
% be larger than window size.
% 
% INPUTS: 
%       x : signal
%       L : Size of the slidding window 
%
% OUTPUTS:
%       y : moving average result
%
  n = length(x);
  if n <= L % Exception & Error Handling
    error("The sliding window must be larger than input signal.")
  end
  y = zeros(1, n);
  for i = 1:n
      if i < L % if the window cannot be centered
          y(i) = mean(x(1:i));
      else % if the window can be centered
          y(i) = mean(x(i-L+1:i));
      end
  end
end


%% References
% [1] https://www.mathworks.com/help/matlab/ref/movmean.html?s_tid=doc_ta
% [2] https://www.mathworks.com/help/dsp/ug/sliding-window-method-and-exponential-weighting-method.html#bvdodb7 
% [3] A. V. Oppenheim, R. W. Schafer, and J. R. Buck, Discrete-Time Signal Processing, 3rd. ed, Prentice Hall Press, 2009.