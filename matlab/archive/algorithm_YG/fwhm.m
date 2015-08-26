function [centroid, width] = fwhm(x, y)
% find FWHM of a distribution 
%   x: index or bins
%   y: value or counts
%   centroid:  x index that gives the maximum y 
if any(y)
    [maxVal, maxIdx] = max(y);
    centroid = maxIdx;
    half = maxVal/2;
    i = 1;
    while sign(y(i)-half) == sign(y(i+1)-half)
        i = i + 1;
    end
    
    xLeft = interp1(y(i:i+1), x(i:i+1), half);
    leftIdx = i;
    i = i + 1;
    while sign(y(i)-half) == sign(y(i+1)-half)
        i = i + 1;
    end
    xRight = interp1(y(i:i+1), x(i:i+1), half);
    rightIdx = i;
    width = xRight - xLeft;
else
    centroid = nan;
    width = inf;
    

end