function [aunc,bunc] = ParametrizeUncertainties(da,bMeas,bTrue, plotflag)
%function [aunc,bunc] = ParametrizeUncertainties(da,bMeas,bTrue, plotflag)
%
% Alpha and beta uncertainties. 10/23/2012.
%
% da and bmeas and btrue should be clean vectors (NaN and Inf removed).
%
% bmeas is used to filter out bmeas=0 as a separate component.
%
% aunc:
%   peakEfficiency: fraction of events in peak around da=0.
%   FWHM: full width at half maximum of peak around da=0.
%   randomEfficiency: fraction of events distributed evenly in da.
%   backscatterEfficiency: fraction of events in peak around abs(da)=180.
%   backscatterFWHM: full width at half maximum of peak around abs(da)=180.
%   RMS: RMS value of alpha error. (added 12/18/12)
%   err.FWHM: 1-sigma uncertainty in FWHM
%   err.peakEfficiency: 1-sigma uncertainty in peak efficiency
%   err. ... etc.
%   x: x values of alpha distribution
%   n: y values of alpha distribution
%   fit: fit object for (x,n) fit
%
% bunc:
%   zeroEfficiency: fraction events at bMeas = 0.
%   FWHM: full width at half maximum of non-zero events, simply measured on the data.
%   position: peak position of non-zero events.
%   peakEfficiency: fraction of events that are non-zero. (1 - zeroEfficiency)
%   RMS: RMS value of beta error. (this was added 12/18/12)
%   err.FWHM: 1-sigma uncertainty in FWHM
%   err.position: 1-sigma uncertainty in non-zero peak position
%   x: x values of beta distribution
%   n: y values of beta distribution
%   fit: fit object for (x,n) fit
%
% plotflag = true;    %testing

%clean alpha
if any(isnan(da)) || any(isinf(da))
    warning('Found NaN''s or Inf''s in data . . removing!')
    da = da(~isnan(da) & ~isinf(da));
end
%clean beta
if length(bMeas) ~= length(bTrue)
    error('bMeas and bTrue must be the same length')
elseif any(isnan(bMeas)) || any(isinf(bMeas)) || any(isnan(bTrue)) || any(isinf(bTrue))
    lg = ~isnan(bMeas) & ~isinf(bMeas) & ~isnan(bTrue) & ~isinf(bTrue);
    bMeas = bMeas(lg);
    bTrue = bTrue(lg);
    warning('Found NaN''s or Inf''s in data . . removing!')
end

Na = length(da);
Nb = length(bMeas);
if length(bMeas) ~= Na
    warning('Vectors da, bMeas are not the same length')
end

%double check range of da.
da(da>180) = da(da>180) - 360;
da(da<-180) = da(da<-180) + 360;

%never mind the sign of da
da = abs(da);

%% Alpha.
alphaResolution = min(100*180/Na, 15);  %histogram resolution scales with statistics
numBins = 180/alphaResolution;
[na,xa] = hist(da,  linspace(alphaResolution/2,180-alphaResolution/2,numBins));

%first, measure alpha FWHM by hand to estimate sigma for fit
alphaMax = max(na);
alphaHalfWidthRightInd = find(na>alphaMax/2, 1, 'last');      %index on the left of the crossover
%linear interpolation
%   left side is at edge of histogram
alphaHalfWidthLeftValue = xa(1);
if alphaHalfWidthRightInd < length(na)
    alphaHalfWidthRightValue = ...
        interp1(na(alphaHalfWidthRightInd + [0,1]), ...
        xa(alphaHalfWidthRightInd + [0,1]), ...
        alphaMax/2);
else
    %up against right side of the distribution. just take index value
    alphaHalfWidthRightValue = xa(alphaHalfWidthRightInd);
end
alphaFWHM_measured = alphaHalfWidthRightValue - alphaHalfWidthLeftValue;
%if half max is in the background statistical noise, this will not work well..
% sigmaMin = 0.6*(alphaFWHM_measured / 2.355);
sigmaMin = alphaResolution;
% sigmaMax = 1.5*(alphaFWHM_measured / 2.355);

%real quantities come from the fit
%alpha fit equation
alphaFitEquation = 'c + a1*exp(-(x.^2./(2*s1^2))) + a2*exp(-(abs(x)-180).^2./(2*s2^2))';
%Variables: a1, a2, c, s1, s2
alphaFitOptions = fitoptions(alphaFitEquation, 'Lower', [0, 0, min(na), sigmaMin, 30], ...
    'Upper', [1.1*max(na), 1.1*max(na), 1.1*mean(na), 100, 100], ...
    'StartPoint', [0.9*max(na(xa<90)), 0.9*max(na(xa>90)), median(na), 40, 50]);
%perform fit
alphaFit = fit(xa',na',alphaFitEquation,alphaFitOptions);

%base values
aunc.peakEfficiency = alphaFit.a1/2 * alphaFit.s1 * sqrt(2*pi) / alphaResolution / Na;     %extra /2 accounts for da = abs(da)
aunc.FWHM = 2.355 * alphaFit.s1;
aunc.backscatterEfficiency = alphaFit.a2/2 * alphaFit.s2 * sqrt(2*pi) / alphaResolution / Na; %extra /2 accounts for da = abs(da)
aunc.randomEfficiency = alphaFit.c * 180 / alphaResolution / Na;
aunc.backscatterFWHM = 2.355 * alphaFit.s2;
%uncertainties go in field "err"
alphaUncertainties = confint(alphaFit,0.68);    %2x5 array of [lower; upper]
aunc.err.peakEfficiency = [aunc.peakEfficiency - (alphaUncertainties(1,1)/2 * alphaUncertainties(1,4) *sqrt(2*pi)/alphaResolution/Na), ...
    (alphaUncertainties(2,1)/2 * alphaUncertainties(2,4) *sqrt(2*pi)/alphaResolution/Na) - aunc.peakEfficiency];
aunc.err.FWHM = [aunc.FWHM - 2.355*alphaUncertainties(1,4), ...
    2.355*alphaUncertainties(2,4) - aunc.FWHM];
aunc.err.backscatterEfficiency = [aunc.backscatterEfficiency - (alphaUncertainties(1,2)/2 * alphaUncertainties(1,5) *sqrt(2*pi)/alphaResolution/Na), ...
    (alphaUncertainties(2,2)/2 * alphaUncertainties(2,5) *sqrt(2*pi)/alphaResolution/Na) - aunc.backscatterEfficiency];
aunc.err.backscatterFWHM = [aunc.backscatterFWHM - 2.355*alphaUncertainties(1,5), ...
    2.355*alphaUncertainties(2,5) - aunc.backscatterFWHM];
aunc.err.randomEfficiency = [aunc.randomEfficiency - alphaUncertainties(1,3)*180/Na, ...
    alphaUncertainties(2,3)*180/alphaResolution/Na];
%save distribution
aunc.x = xa;
aunc.n = na;
%save fit
aunc.fit = alphaFit;
%alpha RMS (independent of other calculations)
aunc.RMS = sqrt(mean(da.^2));

%% Beta zeros.
%first, get bmeas = 0 component.
lgZero = bMeas==0;
bunc.zeroEfficiency = sum(lgZero)/Nb;
bunc.err.zeroEfficiency(1:2) = ones(1,2)*sqrt(Nb*bunc.zeroEfficiency*(1-bunc.zeroEfficiency)) / Nb;
%beta RMS (independent of other calculations)
db = bMeas - abs(bTrue);    %RMS uses both zeros and non-zeros.
bunc.RMS = sqrt(mean(db.^2));
%now, look at remaining data.
db = db(~lgZero);
N2 = sum(~lgZero);
bunc.peakEfficiency = N2/Nb;
bunc.err.peakEfficiency = bunc.err.zeroEfficiency;  %binomial

%do we have enough left?
if N2 < 50
    bunc.FWHM = nan;
    bunc.position = nan;
    bunc.err.FWHM = nan(1,2);
    bunc.err.position = nan(1,2);
    return
end

%% Beta nonzeros.
betaResolution = min(15,15*180/N2);  %need to be able to measure FWHM without too much noise
[nb,xb] = hist(db, -90+betaResolution/2 : betaResolution  : 90-betaResolution/2);

%first, measure beta FWHM by hand to estimate sigma for fit
betaMax = max(na);
betaHalfWidthLeftInd =  find(nb>betaMax/2, 1,'first')-1;    %index on the left of the crossover
betaHalfWidthRightInd = find(nb>betaMax/2, 1, 'last');      %index on the left of the crossover
%linear interpolation
if betaHalfWidthLeftInd > 0
    betaHalfWidthLeftValue = ...
        interp1(nb(betaHalfWidthLeftInd + [0,1]), ...
        xb(betaHalfWidthLeftInd + [0,1]), ...
        betaMax/2);
else
    betaHalfWidthLeftValue = xa(1);
end
if betaHalfWidthRightInd < length(nb)
    betaHalfWidthRightValue = ...
        interp1(nb(betaHalfWidthRightInd + [0,1]), ...
        xb(betaHalfWidthRightInd + [0,1]), ...
        betaMax/2);
else
    %up against right side of the distribution. just take index value
    betaHalfWidthRightValue = xb(betaHalfWidthRightInd);
end
betaFWHM_measured = betaHalfWidthRightValue - betaHalfWidthLeftValue;
%if half max is in the background statistical noise, this will not work well..
% sigmaMin = 0.6*(betaFWHM_measured / 2.355);
sigmaMin = betaResolution;
% sigmaMax = 1.5*(alphaFWHM_measured / 2.355);

%real quantities come from the fit
%beta fit equation
betaFitEquation = 'a1*exp(-((x-b1).^2./(2*s1^2)))';  %single gaussian
%Variables: a1, b1, s1
betaFitOptions = fitoptions(betaFitEquation, 'Lower', [0.25*max(nb), min(db), sigmaMin], ...
    'Upper', [1.2*max(nb), max(db), 40], ...
    'StartPoint', [0.9*max(nb), mean(db), 15]);
%perform fit
betaFit = fit(xb',nb',betaFitEquation,betaFitOptions);

%base values
bunc.FWHM = 2.355 * betaFit.s1;
bunc.position = betaFit.b1;
%uncertainties go in field "err"
betaUncertainties = confint(betaFit,0.68);    %2x3 array of [lower; upper]
bunc.err.FWHM = [bunc.FWHM - 2.355*betaUncertainties(1,3), ...
    2.355*betaUncertainties(2,3) - bunc.FWHM];
bunc.err.position = [bunc.position - betaUncertainties(1,2), ...
    betaUncertainties(2,2) - bunc.position];
%save distribution
bunc.x = xb;
bunc.n = nb;
%save fit
bunc.fit = betaFit;

%% testing
if plotflag
    figure; 
    pos = get(gcf,'position');
    newpos = [pos(1:2),pos(3)*2,pos(4)];
    set(gcf,'position',newpos);
    subplot(1,2,1);
    hold on;
    plot(xa,na,'k');
    plot(alphaFit);
    xlim([0,180]);
    subplot(1,2,2);
    hold on;
    plot(xb,nb,'b');
    plot(betaFit);
    xlim([-90,90]);
    drawnow;
    disp([]);   %debug breakpoint
end