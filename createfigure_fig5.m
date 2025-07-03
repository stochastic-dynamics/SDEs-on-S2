function createfigure_fig5(X1, YMatrix1)
%CREATEFIGURE(X1, YMatrix1)
%  X1:  vector of plot x data
%  YMATRIX1:  matrix of plot y data
%
% Create figure
figure('OuterPosition',[31 12 1362 1039]);
%
% Create axes
axes1 = axes('Position',[0.13 0.12565445026178 0.775 0.799345549738222]);
hold(axes1,'on');
%
% Create multiple line objects using matrix input to plot
plot1 = plot(X1,YMatrix1,'Marker','o','LineWidth',2);
set(plot1(1),'DisplayName','g-Weak 1.0','Color',[1 0 0]);
set(plot1(2),'DisplayName','g-Weak 2.0','Color',[0 1 0]);
set(plot1(3),'DisplayName','g-Weak 3.0','Color',[0 0 1]);
%
% Create ylabel
ylabel('$\|\bf{E}[z_T]-E[\bar{z}_T]\|$','Interpreter','latex');
%
% Create xlabel
xlabel('\Deltat');
%
% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[(-3) (0)]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[(-5) (-1)]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontName','Times','FontSize',24,'FontWeight','bold','XMinorGrid',...
    'on','XTick',[-3 -2 -1 0],'XTickLabel',...
    {'10^{-3}','10^{-2}','10^{-1}','10^{0}'},'YMinorGrid','on','YTick',...
    [-5 -4 -3 -2 -1],'YTickLabel',...
    {'10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}'});
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.140941010205558 0.802468741815967 0.141901928549536 0.109936571934006]);
%
