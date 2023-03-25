%% [TIP7200 - Processamento Digital de Sinais] - Moving Average
% Authour: Lucas Abdalah
% ----------------------------
% 
% moving_average.m
% 2023/03/24 - v1
%

%% Prepare Ambient 
clearvars;
close all; 


%% Signal Input
x = [1, 2, -1, 2, 1, 1, -2, 1];
L = 5;
shifted = (1:length(x)) + floor(L/2);
y = ma(x, L);
y_matlab = movmean(x, L);


%% Visualize Moving Average
h1 = figure();
subplot(2,1,1); % Signal vs. implemented MA  vs. Matlab MA
stem(x,...
  'Color', 'black',...        
  'LineStyle', ':',...
  'LineWidth', 1.5);
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
ylabel("x[n]")
xlim([-1, 11])
legend("x", "y", "y\_matlab", 'Location', 'northeastoutside','Orientation', 'Horizontal','Position', [0.5 0.47 0.0 1], 'Units','normalized');
legend boxoff

subplot(2,1,2); % Signal vs. shifted implemented MA  vs. Matlab MA
stem(x,...
  'Color', 'black',...        
  'LineStyle', ':',...
  'LineWidth', 1.5);
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
xlim([-1, 11])
xlabel("n")
ylabel("x[n]")
annotation('textarrow',[0.48 0.42], [0.18 0.32],'String',"Sample shift [+2]")
annotation('rectangle',[.25 .30 .1 .13],'Color','red')
annotation('rectangle',[.76 .24 .1 .13],'Color','blue')


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
% [1] https://www.mathworks.com/help/dsp/ug/sliding-window-method-and-exponential-weighting-method.html#bvdodb7 