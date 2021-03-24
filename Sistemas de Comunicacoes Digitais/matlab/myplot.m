function XNormalized = myplot(x,y)
    % Computes z-score data normalization (and standardization) of input signal(s)
    %
    % SYNTAX: XNormalized = center_scale(X);
    %
    %         XNormalized : input signals with zero mean and std one
    %
    %         X  : input signals (one signal per row, one sample per column).
    %
    % HISTORY:
    %
    % 2019/04/30: - created by Lucas Abdalah.


%% PCA Plot - 2D plan
figure
s1 = scatter(score(:,1),score(:,2),'filled');
set(gca,'Color','none')
s1.MarkerFaceColor = 'b';
text(score(:,1),score(:,2),labels,...
     'VerticalAlignment','bottom',...
     'HorizontalAlignment','center')
title('2D - Principals Components Plan')
xlabel('pc1')
ylabel('pc2')
grid on