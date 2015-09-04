function printAll(figureHandle,fileName)
% function printAll(figureHandle,fileName)
%
% Save figure represented by handle figureHandle to files in
%   .eps, .png and .fig formats.
%
% figureHandle can be gcf for the current figure.
%
% fileName does not need to include the .fig suffix.

%do not duplicate the file suffix if it is already given
if strcmpi(fileName(end-3:end), '.fig') || ...
        strcmpi(fileName(end-3:end), '.png') || ...
        strcmpi(fileName(end-3:end), '.eps')
    fileName = fileName(1:end-4);
end

%eps
print(figureHandle,'-painters','-depsc2',[fileName,'.eps']);

%png
print(figureHandle,'-painters','-dpng',[fileName,'.png']);

%fig
saveas(figureHandle,fileName,'fig');