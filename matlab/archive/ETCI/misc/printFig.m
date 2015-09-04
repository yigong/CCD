function printFig(figureHandle,fileName)
% function printFig(figureHandle,fileName)
%
% Save figure represented by handle figureHandle to a file, fileName.fig.
%
% figureHandle can be gcf for the current figure.
%
% fileName does not need to include the .fig suffix.

%do not duplicate the file suffix if it is already given
if strcmpi(fileName(end-3:end), '.fig')
    saveas(figureHandle,fileName);
else
    saveas(figureHandle,fileName,'fig');
end
