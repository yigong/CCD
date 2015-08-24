function [Centroid, Std] = noiseFitSeg(E, lowerLimit, upperLimit)
% function [Centroid, Std] = noiseFitSeg(E, lowerLimit, upperLimit)
%
% Fit data E in range given by lowerLimit and upperLimit with a weighted
% 'gauss1' fit. It outputs the centroid and standard deviation of the
% gaussian fit.
%
% Input:
% lowerLimit: lowerlimit of the range to fit, in ADC unit.
% upperLimit: upperlimit of the range to fit, in ADC unit.
% E: data to fit, in ADC unit
% 
% Output:
% Centroid: in ADC unit
% Std: in ADC unit

EtoFit = E(E > lowerLimit & E < upperLimit);
nBin_Rice = ricerule(EtoFit);
binWidth_temp = (upperLimit - lowerLimit)/nBin_Rice;
binWidth = floor(binWidth_temp/0.5)*0.5;
[y,x] = hist(EtoFit,lowerLimit:binWidth:upperLimit);
weights = 1./(y+0.25);
fitOpt = fitoptions('gauss1','Weights', weights,...
    'StartPoint',[max(y),median(x),(upperLimit-lowerLimit)/3]);
noiseFit = fit(x', y', 'gauss1', fitOpt);
Centroid = noiseFit.b1;
Std = noiseFit.c1/sqrt(2);

end