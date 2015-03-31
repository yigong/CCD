function FixAxesMargins(h_axes)
%function FixAxesMargins(h_axes)
%
% Add vertical space inside figure, above and below axes, in order to 
% get x axis label to show.
%
% h_axes: handle of axes

axes_dposition = [0.05,0.1,-0.1,-0.15];     %units = normalized

% Note: position = [left bottom width height]

h_figure = get(h_axes,'Parent');

lastAxesUnits = get(h_axes,'Units');
set(h_axes,'Units','Normalized');
set(h_axes,'ActivePositionProperty','OuterPosition')

axes_position = get(h_axes,'Position');

set(h_axes,'Position',[axes_position + axes_dposition]);

lastFigureUnits = get(h_figure,'Units');
set(h_figure,'Units','Pixels');

figure_position = get(h_figure,'Position');
new_figure_position = ...
    [figure_position(1) - axes_dposition(1)*figure_position(3), ...
    figure_position(2) - axes_dposition(2)*figure_position(4), ...
    figure_position(3) + abs(axes_dposition(3))*figure_position(3), ...
    figure_position(4) + abs(axes_dposition(4))*figure_position(4)];
set(h_figure,'Position',new_figure_position);

%return to previous units
set(h_axes,'Units',lastAxesUnits);
set(h_figure,'Units',lastFigureUnits);
%old way doesn't accomodate resizing of figure.
%{
dheight_figure = 100;
dbottom_axes = 50;



h_figure = get(h_axes,'Parent');

set(h_axes,'Units','pixels');   %remove normalization

axes_position = get(h_axes,'Position');
figure_position = get(h_figure,'Position');

%make figure box taller
set(h_figure,'Position',[figure_position(1:3), figure_position(4) + dheight_figure]);

%make axes box higher
set(h_axes,'Position',[axes_position(1), axes_position(2) + dbottom_axes, axes_position(3:4)]);
%}