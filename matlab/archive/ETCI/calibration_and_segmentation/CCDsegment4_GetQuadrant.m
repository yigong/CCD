function [xOffset, yOffset, imageQuadrant] = CCDsegment4_GetQuadrant(dim1, dim2, imageFull, quadrantIndex)
%function [xOffset, yOffset, imageQuadrant] = CCDsegment4_GetQuadrant(dim1, dim2, imageFull, quadrantIndex)
%
% Operated by CCDsegment4.

if dim1/2 ~= round(dim1/2) || dim2/2 ~= round(dim2/2)
    error('Image dimensions must be even in order to divide into quadrants.')
end

switch quadrantIndex
    %Quadrants are numbered in clockwise order from upper left
    case 1
        %bottom
        xQuadrant = 1:dim1/2;
        
    case 2
        %top
        xQuadrant = dim1/2+1:dim1;
        
end

%grab the image of one quadrant
imageQuadrant = imageFull(xQuadrant,:);

xOffset = xQuadrant(1);
yOffset = 1;
