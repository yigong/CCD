function createfigure1(X1, YMatrix1)
%CREATEFIGURE1(X1, YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 23-Jul-2015 15:54:27

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'FontSize',18);
box(axes1,'on');
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'LineWidth',2,'Parent',axes1);
set(plot1(1),'DisplayName','100nm','Color',[0 0 1]);
set(plot1(2),'DisplayName','100nm','LineStyle','--','Color',[0 0 1]);
set(plot1(3),'DisplayName','10um','Color',[1 0 0]);
set(plot1(4),'DisplayName','10um','LineStyle','--','Color',[1 0 0]);
set(plot1(5),'DisplayName','1um','Color',[0 1 0]);
set(plot1(6),'DisplayName','1um','LineStyle','--','Color',[0 1 0]);
set(plot1(7),'DisplayName','300nm','Color',[0 1 1]);
set(plot1(8),'DisplayName','300nm','LineStyle','--','Color',[0 1 1]);
set(plot1(9),'DisplayName','3um','Color',[1 0 1]);
set(plot1(10),'DisplayName','3um','LineStyle','--','Color',[1 0 1]);
set(plot1(11),'DisplayName','50um','Color',[0 0 0]);
set(plot1(12),'DisplayName','50um','LineStyle','--','Color',[0 0 0]);

% Create xlabel
xlabel('mm','FontSize',18);

% Create ylabel
ylabel('counts','FontSize',18);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'FontSize',18);

