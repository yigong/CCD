function threshold = CCDsegment4_GetPreCalThreshold(imagePortion,thresholdSigma)
%function threshold = CCDsegment4_GetPreCalThreshold(imagePortion,thresholdSigma)
%
% Operated by CCDsegment4.

%the threshold is given relative to the noise (in sigma) of the black level.
%we histogram all the pixels, which should be completely dominated by
%  black level instead of tracks, and fit this histogram to get the 
%  black level noise width. That tells us the numerical value of the 
%  threshold.
imageMedian = median(imagePortion(:)); %median should be 0 or very close to 0
imageStDev = std(imagePortion(:));    %standard deviation
%make a histogram of pixel values in the range of +10/-6 st. dev.'s
[n,x] = hist(imagePortion(:), (imageMedian-6*imageStDev):(imageStDev/10):(imageMedian+10*imageStDev));
%perform a fit
blackLevelFitOptions = fitoptions('Method','LinearLeastSquares',...
    'Upper',[1.5*max(n),imageMedian+1,sqrt(2)*imageStDev*1.2],...
    'Lower',[0.8*max(n),imageMedian-1,sqrt(2)*imageStDev*0.8]);
blackLevelFit = fit(x',n','gauss1',blackLevelFitOptions);

%threshold is determined by the width of the gaussian fit.
threshold = imageMedian + thresholdSigma * blackLevelFit.c1 / sqrt(2);