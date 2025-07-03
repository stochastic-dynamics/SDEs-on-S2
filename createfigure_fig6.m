function createfigure_fig6(Y1, Y2, Y3)
%CREATEFIGURE(Y1, Y2, Y3)
%  Y1:  vector of plot y data
%  Y2:  vector of plot y data
%  Y3:  vector of plot y data
%
% Create figure
figure;
%
% Create subplot
subplot1 = subplot(3,1,1);
hold(subplot1,'on');
%
% Create plot
plot(Y1,'LineWidth',2);
%
% Create ylabel
ylabel('q_1');
%
xlim(subplot1,[0 1000]);
box(subplot1,'on');
grid(subplot1,'on');
hold(subplot1,'off');
% Set the remaining axes properties
set(subplot1,'FontName','Times New Roman','FontSize',24,'FontWeight','bold',...
    'XTick',[0 200 400 600 800 1000],'XTickLabel',...
    {'0','200','400','600','800','1000'});
% Create subplot
subplot2 = subplot(3,1,2);
hold(subplot2,'on');
%
% Create plot
plot(Y2,'LineWidth',2);
%
% Create ylabel
ylabel('q_2');
%
xlim(subplot2,[0 1000]);
box(subplot2,'on');
grid(subplot2,'on');
hold(subplot2,'off');
% Set the remaining axes properties
set(subplot2,'FontName','Times New Roman','FontSize',24,'FontWeight','bold',...
    'XTick',[0 200 400 600 800 1000],'XTickLabel',...
    {'0','200','400','600','800','1000'},'YTick',[-0.2 0.2 0.6]);
% Create subplot
subplot3 = subplot(3,1,3);
hold(subplot3,'on');
%
% Create plot
plot(Y3,'LineWidth',2);
%
% Create ylabel
ylabel('q_3');
%
% Create xlabel
xlabel('Samples');
%
xlim(subplot3,[0 1000]);
box(subplot3,'on');
grid(subplot3,'on');
hold(subplot3,'off');
% Set the remaining axes properties
set(subplot3,'FontName','Times New Roman','FontSize',24,'FontWeight','bold',...
    'XTick',[0 200 400 600 800 1000],'XTickLabel',...
    {'0','200','400','600','800','1000'});
%
% End