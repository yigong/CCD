function x = SampleGauss2(varargin)
%x = SampleGauss2
%x = SampleGauss2(N)
%x = SampleGauss2(N,M)
%x = SampleGauss2([N,M])
%
%Sample a gaussian distribution of sigma=1 and mean=0.
%
% N = length of vector of values to return
% N,M = row,col of array of values to return
% [N,M] = row,col of array of values to return (e.g. size(A))
% If no arguments are specified, returns a single value.

if nargin==2
    %SampleGauss2(N,M)
    r = rand(varargin{1},varargin{2});
elseif nargin==1 && length(varargin{1})==1
    %SampleGauss2(N)
    r = rand(1,varargin{1});
elseif nargin==1 && length(varargin{1})==2
    %SampleGauss2([N,M])
    r = rand(varargin{1}(1),varargin{1}(2));
elseif nargin==0
    %SampleGauss2
    r = rand;
else
    x = NaN;
    return
end

%input to erfinv should be on [-1,1]
%sqrt(2) adjusts sigma
x = sqrt(2)*erfinv(r.*2-1);