function f = fit_copy(varargin)
%function f = fit_copy(xData, yData)
%
% Take gaussian-like data, and measure a FWHM/1.665 value.
% This is similar to what a gaussian fit would result in, but without any
%   actual fitting.
%
 
% %measure centroid
% xCentroid = sum(xData.*yData) / sum(yData);
% 
% %find local maximum closest to centroid?
% %...
% 

xData = varargin{1};
yData = varargin{2};
%ignore the rest

%find global maximum
[valMax, indMax] = max(yData);

%find sides of FHWM
fLeft =               find(yData(2:indMax)     > valMax/2 & yData(1:indMax-1)   < valMax/2, 1, 'last');   %left index of crossover point
fRight = indMax - 1 + find(yData(indMax:end-1) > valMax/2 & yData(indMax+1:end) < valMax/2, 1, 'first');  %left index of crossover point
% 
% %ignore any crossovers on the wrong side
% fLeft =   fLeft(fLeft < indMax);
% fRight = fRight(fRight > indMax);
% 
% %just take the crossover closest to the middle
% fLeft = fLeft(end);
% fRight = fRight(1);

%problematic cases
if isempty(fLeft) || isempty(fRight)
    %uhhh . . need to take some value or else alpha,beta calculation will fail completely
    f.c1 = 0;
else
    %linear interpolation to get the half-max point
    pos1 = xData(fRight) + (valMax/2 - yData(fRight)) / (yData(fRight+1) - yData(fRight)) * (xData(fRight+1) - xData(fRight));
    pos2 = xData(fLeft)  + (valMax/2 - yData(fLeft))  / (yData(fLeft +1) - yData(fLeft))  * (xData(fLeft +1) - xData(fLeft));
    f.c1 = (pos1 - pos2) / 2.355 * sqrt(2);    %matlab definition of c1 is sqrt(2)*sigma
    %   units are image pixels
end