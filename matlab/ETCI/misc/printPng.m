function printPng(figureHandle,fileName)
% function printPng(figureHandle,fileName)
%
% Save figure represented by handle figureHandle to a file, fileName.png.
%
% figureHandle can be gcf for the current figure.
%
% fileName does not need to include the .png suffix.

%do not duplicate the file suffix if it is already given
if strcmpi(fileName(end-3:end), '.png')
    print(figureHandle,'-painters','-dpng',fileName);
else
    print(figureHandle,'-painters','-dpng',[fileName,'.png']);
end
