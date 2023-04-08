%% [TIP7200 - Processamento Digital de Sinais] - Convolution Sum
% Authour: Lucas Abdalah
% ----------------------------
% 
% ex3_conv_sum.m
% 2023/03/29 - v1
% 

%% Prepare Ambient 
clearvars;
close all;
% pause on

%% Synthetic signal Input
x = [1, 2, -1, 2, 1, 1, -2, 1];
n = 1:length(x);
h = [1 0.5 0 -0.25]; % LTI system h[n] 
y_matlab = conv(x, h, "same");
y = conv_sum(x, h);

figure(1);
stem(n, x, 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5);
hold on;
plot(n, y(3:10), ...
  'Color', 'red',...        
  'LineStyle', '-',...
  'LineWidth', 1.5,...
  'Marker', 'v',...
  'MarkerFaceColor', 'red',...
  'MarkerSize', 5);
plot(n, y_matlab,...
  'Color', 'blue',...        
  'LineStyle', '--',...
  'LineWidth', 1.5,...
  'Marker', 'square',...
  'MarkerFaceColor', 'blue',...
  'MarkerSize', 5);
hold off;
legend("x", "y", "y\_matlab", 'Location', 'northeastoutside','Orientation', 'Horizontal','Position', [0.5 0.47 0.0 1], 'Units','normalized');
legend boxoff
xlabel("n");
ylabel("Amplitude");
grid on;


% %% Audio Signal Input 1
pause()
[x_fala, Fs_fala] = audioread("fala_sino.wav");
N_fala = length(x_fala);
t_fala = linspace(0, N_fala/Fs_fala, N_fala);
y_fala = conv_sum(x_fala, h);
y_fala = transpose(y_fala(3:end-1));
window = 1:250;

filename = ("../audio/ex3_y_fala.wav");
audiowrite(filename, y_fala, Fs_fala);

figure(2);
plot(t_fala(window), x_fala(window), 'Color', 'black', 'LineStyle', '-', 'LineWidth', 1.5);
hold on;
plot(t_fala(window), y_fala(window), 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.25);
hold off;
legend("x", "y", 'Location', 'northeastoutside','Orientation', 'Horizontal','Position', [0.5 0.47 0.0 1], 'Units','normalized');
legend boxoff
xlabel("Time (s)")
ylabel("Amplitude")
grid on;
axis tight


%% Audio Signal Input 2
[x_cantina, Fs_cantina] = audioread("cantinaband.wav");
x_cantina = x_cantina.';
N_cantina = length(x_cantina);
t_cantina = linspace(0, N_cantina/Fs_cantina, N_cantina);
y_cantina = conv_sum(x_cantina, h);
y_cantina = transpose(y_cantina(3:end-1));
window = 1:250;

filename = ("../audio/ex3_y_cantina.wav");
audiowrite(filename, y_cantina, Fs_cantina);

figure(3);
plot(t_cantina(window), x_cantina(window), 'Color', 'black', 'LineStyle', '-', 'LineWidth', 1.5);
hold on;
plot(t_cantina(window), y_cantina(window), 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
hold off;
legend("x", "y", 'Location', 'northeastoutside','Orientation', 'Horizontal','Position', [0.5 0.47 0.0 1], 'Units','normalized');
legend boxoff
xlabel("Time (s)")
ylabel("Amplitude")
grid on;
axis tight


%% Save results
% savefig_tight(figure(1), '../figures/ex3_synthetic', 'both');
% savefig_tight(figure(2), '../figures/ex3_fala_sino', 'both');
% savefig_tight(figure(3), '../figures/ex3_cantina', 'both');


%% Local Functions
function y = conv_sum(x, h)

  % Check if x and h are row vectors
  x = checkTranspose(x);
  h = checkTranspose(h);
  
  % Compute the length of the input and filter sequences
  Nx = length(x);
  Nh = length(h);

  y = zeros(1, Nx + Nh - 1);

  for n = 1:Nx
      for k = 1:Nh
          y(n+k-1) = y(n+k-1) + x(n)*h(k);
      end
  end
  
end


function y = checkTranspose(x)
  % Check the size of the input vector
  [m, n] = size(x);

  % If the vector is column-wise, transpose it
  if n == 1 && m > 1
      y = x.';
  else
      y = x;
  end
end

%% References
% [1] https://www.mathworks.com/help/matlab/ref/conv.html
% [2] A. V. Oppenheim, R. W. Schafer, and J. R. Buck, Discrete-Time Signal Processing, 3rd. ed, Prentice Hall Press, 2009.