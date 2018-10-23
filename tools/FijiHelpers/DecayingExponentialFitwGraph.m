function cf_ = DecayingExponentialFitwGraph(time_vect,bead1)
%CREATEFIT Create plot of data sets and fits
%   CREATEFIT(TIME_VECT,BEAD1)
%   Creates a plot, similar to the plot in the main Curve Fitting Tool,
%   using the data that you provide as input.  You can
%   use this function with the same data you used with CFTOOL
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of data sets:  2
%   Number of fits:  1

% Data from data set "bead1 vs. time_vect":
%     X = time_vect:
%     Y = bead1:
%     Unweighted

% Data from data set "bead1 vs. time_vect (2 )":
%     X = time_vect:
%     Y = bead1:
%     Unweighted

% Auto-generated by MATLAB on 03-Jun-2015 16:19:21

% Set up figure to receive data sets and fits
f_ = clf;
% figure(f_);
% f_ = f0;
set(f_,'Units','Pixels','Position',[1037 398 672 475]);
% Line handles and text for the legend.
legh_ = [];
legt_ = {};
% Limits of the x-axis.
xlim_ = [Inf -Inf];
% Axes for the plot.
ax_ = axes;
set(ax_,'Units','normalized','OuterPosition',[0 0 1 1]);
set(ax_,'Box','on');
axes(ax_);
hold on;

% --- Plot data that was originally in data set "bead1 vs. time_vect"
time_vect = time_vect(:);
bead1 = bead1(:);
h_ = line(time_vect,bead1,'Parent',ax_,'Color',[0.333333 0 0.666667],...
    'LineStyle','none', 'LineWidth',1,...
    'Marker','.', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(time_vect));
xlim_(2) = max(xlim_(2),max(time_vect));


% --- Plot data that was originally in data set "bead1 vs. time_vect (2 )"
h_ = line(time_vect,bead1,'Parent',ax_,'Color',[0.333333 0.666667 0],...
    'LineStyle','none', 'LineWidth',1,...
    'Marker','.', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(time_vect));
xlim_(2) = max(xlim_(2),max(time_vect));
legh_(end+1) = h_;
legt_{end+1} = 'bead1 vs. time_vect (2 )';

% Nudge axis limits beyond data limits
if all(isfinite(xlim_))
    xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
    set(ax_,'XLim',xlim_)
else
    set(ax_, 'XLim',[-0.11000000000000000056, 11.109999999999999432]);
end

% --- Create fit "fit 1"
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[-10 -Inf]);
ok_ = isfinite(time_vect) & isfinite(bead1);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs',...
        'Ignoring NaNs and Infs in data.' );
end
st_ = [bead1(end) time_vect(2) ];
set(fo_,'Startpoint',st_);
ft_ = fittype('a*(-exp(-x/tau)+1)',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a', 'tau'});

% Fit this model using new data
cf_ = fit(time_vect(ok_),bead1(ok_),ft_,fo_);
% Alternatively uncomment the following lines to use coefficients from the
% original fit. You can use this choice to plot the original fit against new
% data.
%    cv_ = { 1.0747077815443111248, 2.7136877198306561532};
%    cf_ = cfit(ft_,cv_{:});

% Plot this fit
h_ = plot(cf_,'fit',0.95);
set(h_(1),'Color',[1 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
% Turn off legend created by plot method.
legend off;
% Store line handle and fit name for legend.
legh_(end+1) = h_(1);
legt_{end+1} = 'fit 1';

% --- Finished fitting and plotting data. Clean up.
hold off;
% Display legend
leginfo_ = {'Orientation', 'vertical'};
h_ = legend(ax_,legh_,legt_,leginfo_{:});
set(h_,'Units','normalized');
t_ = get(h_,'Position');
t_(1:2) = [0.498884,0.401579];
set(h_,'Interpreter','none','Position',t_);
% Remove labels from x- and y-axes.
xlabel(ax_,'');
ylabel(ax_,'');
