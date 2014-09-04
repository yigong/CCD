function [aunc,bunc] = CalculateUncertainties(da,db,conflevel)
%function [aunc,bunc] = CalculateUncertainties(da,db,conflevel)
%
% Alpha and beta uncertainties. 6/9/2010.
% conflevel = half-angle confidence level. E.g. 0.68 for 68%.
%
% da and db should be clean vectors (NaN and Inf removed).
%
% aunc, bunc are each 3-element vectors.
% aunc(1) is the measured half-angle uncertainty.
% aunc(2) is the lower value of the error bar (corresponding to Nlow = N0 - sqrt(variance_binomial)).
% aunc(3) is the upper value of the error bar ( " ).
%
% similarly for bunc(1:3).

res = 0.01;     %histogram resolution, in degrees. for both alpha and beta.

N = length(da);
if N ~= length(db)
    error('Vectors da, db must be the same length')
end

%check for events
if N < 5
    warning('Too few events to calculate uncertainties...')
    aunc = nan(1,3);
    bunc = nan(1,3);
    return
end

N0 = N * conflevel;

Nsig = sqrt(N * conflevel * (1-conflevel)); %variance of the binomial distribution

%values to use for the error bar
Nlow = N0 - Nsig;
Nhigh = N0 + Nsig;

%double check range of da.
da(da>180) = da(da>180) - 360;
da(da<-180) = da(da<-180) + 360;

%never mind the sign of da, db.
da = abs(da);
db = abs(db);

[na,xa] = hist(da,0:res:180);
[nb,xb] = hist(db,0:res:90);

na = cumsum(na);
nb = cumsum(nb);

%half-angle uncertainty, and location of error bars.
aunc(1) = xa(find(na>N0,1));
aunc(2) = xa(find(na>Nlow,1));
aunc(3) = xa(find(na>Nhigh,1));

bunc(1) = xb(find(nb>N0,1));
bunc(2) = xb(find(nb>Nlow,1));
bunc(3) = xb(find(nb>Nhigh,1));

%convert to format used by errorbar function
%errorbar(X,Y,L,U,...)
% L is the length of the errorbar below the point; U is the length above

aunc(2) = aunc(1) - aunc(2);
aunc(3) = aunc(3) - aunc(1);

bunc(2) = bunc(1) - bunc(2);
bunc(3) = bunc(3) - bunc(1);