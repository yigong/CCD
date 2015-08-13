function [centroid, width] = fwhm(x, y)
% find FWHM of a distribution 
%   x: index or bins
%   y: value or counts
%   centroid:  x index that gives the maximum y 
    [maxVal, maxIdx] = max(y);
    centroid = maxIdx;
    half = maxVal/2;
    i = 1;
    while sign(y(i)-half) == sign(y(i+1)-half)
        i = i + 1;
    end
    leftIdx = i;
    i = i + 4;
    while sign(y(i)-half) == sign(y(i+1)-half)
        i = i + 1;
    end
    rightIdx = i;
    width = x(rightIdx) - x(leftIdx);
    
    
end

