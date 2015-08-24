function fit110Kcustom(xdata,ydata)
%FIT110KCUSTOM    Create plot of datasets and fits
%   FIT110KCUSTOM(XDATA,YDATA)
%   Creates a plot, similar to the plot in the main curve fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with cftool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  2
%   Number of fits:  1

 
% Data from dataset "ydata vs. xdata":
%    X = xdata:
%    Y = ydata:
%    Unweighted
 
% Data from dataset "ydata vs. xdata (smooth)":
%    X = xdata:
%    Y = ydata:
%    Unweighted
%
% This function was automatically generated on 06-Aug-2008 14:18:31

% Set up figure to receive datasets and fits
f_ = clf;
figure(f_);
set(f_,'Units','Pixels','Position',[576 189 680 484]);
legh_ = []; legt_ = {};   % handles and text for legend
xlim_ = [Inf -Inf];       % limits of x axis
ax_ = axes;
set(ax_,'Units','normalized','OuterPosition',[0 0 1 1]);
set(ax_,'Box','on');
axes(ax_); hold on;

 
% --- Plot data originally in dataset "ydata vs. xdata"
xdata = xdata(:);
ydata = ydata(:);
% This dataset does not appear on the plot
% Add it to the plot by removing the if/end statements that follow
% and by selecting the desired color and marker
if 0
   h_ = line(xdata,ydata,'Color','r','Marker','.','LineStyle','none');
   xlim_(1) = min(xlim_(1),min(xdata));
   xlim_(2) = max(xlim_(2),max(xdata));
   legh_(end+1) = h_;
   legt_{end+1} = 'ydata vs. xdata';
end       % end of "if 0"
 
% --- Plot data originally in dataset "ydata vs. xdata (smooth)"
sm_.y2 = smooth(xdata,ydata,0.25,'loess',0);
h_ = line(xdata,sm_.y2,'Parent',ax_,'Color',[0.333333 0.666667 0],...
     'LineStyle','none', 'LineWidth',1,...
     'Marker','.', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(xdata));
xlim_(2) = max(xlim_(2),max(xdata));
legh_(end+1) = h_;
legt_{end+1} = 'ydata vs. xdata (smooth)';

% Nudge axis limits beyond data limits
if all(isfinite(xlim_))
   xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
   set(ax_,'XLim',xlim_)
else
    set(ax_, 'XLim',[156.15908364951676, 9849.1741527965205]);
end


% --- Create fit "fit 1"
ok_ = isfinite(xdata) & isfinite(sm_.y2);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end
st_ = [0.23693388099561674 0.31556645868109046 0.15062858815493141 ];
ft_ = fittype('a*x/((1+(x/b)^c)^(1/c))',...
     'dependent',{'y'},'independent',{'x'},...
     'coefficients',{'a', 'b', 'c'});

% Fit this model using new data
cf_ = fit(xdata(ok_),sm_.y2(ok_),ft_,'Startpoint',st_);

% Or use coefficients from the original fit:
if 0
   cv_ = { 4449.1120491440479, 2691.1604838621606, 0.76934879612066231};
   cf_ = cfit(ft_,cv_{:});
end

% Plot this fit
h_ = plot(cf_,'fit',0.95);
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[1 0 0],...
     'LineStyle','-', 'LineWidth',2,...
     'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_(1);
legt_{end+1} = 'fit 1';

% Done plotting data and fits.  Now finish up loose ends.
hold off;
leginfo_ = {'Orientation', 'vertical'}; 
h_ = legend(ax_,legh_,legt_,leginfo_{:}); % create and reposition legend
set(h_,'Units','normalized');
t_ = get(h_,'Position');
t_(1:2) = [0.585784,0.455234];
set(h_,'Interpreter','none','Position',t_);
xlabel(ax_,'');               % remove x label
ylabel(ax_,'');               % remove y label
