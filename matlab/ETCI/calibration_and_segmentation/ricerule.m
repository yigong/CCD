function nBins = ricerule(E)
% function nBins = ricerule(E)
%
% Calculate the number of bins to histgram data based on Rice rule.
%
% Input: 
% E: data to histogram.
% 
% Output:
% nBins: number of bins given by Rice rule.

n = length(E);
nBins = ceil(2*n^(1/3));
end