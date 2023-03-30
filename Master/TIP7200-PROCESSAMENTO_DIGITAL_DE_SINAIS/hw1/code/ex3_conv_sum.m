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
% n = 1:length(x);
h = [1 0.5 0 -0.25]; % LTI system h[n] 
y_matlab = conv(x, h, "same");
y = conv_sum(x, h);

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
plot(y_matlab,...
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


%% Audio Signal Input 1
pause()
[x_fala, Fs_fala] = audioread("fala_sino.wav");
N_fala = length(x_fala);
t_fala = linspace(0, N_fala/Fs_fala, N_fala);
y_fala = conv_sum(x_fala, h);

figure(2);
plot(t_fala, x_fala, 'Color', 'black', 'LineStyle', '-', 'LineWidth', 1.5);
hold on;
plot(t_fala, y_fala, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
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

figure(3);
plot(t_cantina, x_cantina, 'Color', 'black', 'LineStyle', '-', 'LineWidth', 1.5);
hold on;
plot(t_cantina, y_cantina, 'Color', 'red', 'LineStyle', '--', 'LineWidth', 0.5);
hold off;
legend("x", "y", 'Location', 'northeastoutside','Orientation', 'Horizontal','Position', [0.5 0.47 0.0 1], 'Units','normalized');
legend boxoff
xlabel("Time (s)")
ylabel("Amplitude")
grid on;
axis tight


%% Save results
% savefig_tight(figure(1), '../figures/ex2_synthetic', 'both');
% savefig_tight(figure(2), '../figures/ex2_fala_sino', 'both');
% savefig_tight(figure(3), '../figures/ex2_cantina', 'both');


%% Local Functions
function y = conv_sum(x, h)
% function y = conv_sum(x, h, option)

  % Check if x and h are row vectors
  x = checkTranspose(x);
  h = checkTranspose(h);
  
  % Compute the length of the input and filter sequences
  Nx = length(x);
  Nh = length(h);
  % if isempty(option)
  %   option = ""
  % switch option
  % case 'full'
    y = zeros(1, Nx + Nh - 1);

    for n = 1:Nx
        for k = 1:Nh
            y(n+k-1) = y(n+k-1) + x(n)*h(k);
        end
    end
  
  % case 'same'

  % % Compute the length of the convolution and its starting index
  % Ny = Nx + Nh - 1;
  % y_shift = floor((Nh - 1) / 2) + 1;
  % x_pad = [zeros(1, y_shift-1), x, zeros(1, Ny-Nx-y_shift+1)]; % Pad the input sequence with zeros to ensure proper alignment

  % y = zeros(1, Ny);

  % for n = 1:Ny
  %     for k = 1:Nh
  %         if n-k+1 > 0 && n-k+1 <= Nx
  %              y(n) = y(n) + x_pad(n-k+1) * h(k); % Compute the convolution sum using a for loop
  %         end
  %     end
  % end

  % y = y(y_shift:y_shift+Nx-1); % Extract the central part of the convolution that is the same size as x

  % otherwise
  %   error("")
  % end

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