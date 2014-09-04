function printEps(figureHandle,fileName)
% function printEps(figureHandle,fileName)
%
% Save figure represented by handle figureHandle to a file, fileName.eps.
%
% figureHandle can be gcf for the current figure.
%
% fileName does not need to include the .eps suffix.

%do not duplicate the file suffix if it is already given
if strcmpi(fileName(end-3:end), '.eps')
    print(figureHandle,'-painters','-depsc2',fileName);
else
    print(figureHandle,'-painters','-depsc2',[fileName,'.eps']);
end