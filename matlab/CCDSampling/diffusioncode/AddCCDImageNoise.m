function img = AddCCDImageNoise(varargin)
%function img = AddCCDImageNoise(img)
%function img = AddCCDImageNoise(img,'blsig',blsig)
%function img = AddCCDImageNoise(img,'blfwhm',blfwhm)
%function img = AddCCDImageNoise(img,'fano',fano)
%function img = AddCCDImageNoise(img,'epcc',epcc)
%
% Add statistical noise to the pixel values of a simulated CCD image.
%
% Inputs:
%   img:    the ideal CCD image, with pixel values in keV
% Optional input arguments: [default value]
%   blsig:  the 1-sigma of the black level (fixed) image noise, in keV [0.01874]
%   blfwhm: the FWHM of the black level (fixed) image noise, in keV [2.355*0.01874]
%   fano:   the Fano factor to use for the carrier statistics [0.14]
%   epcc:   the average energy to create a charge carrier (electron-hole pair), in eV [3.7]
%
% Output:
%   img:    the noisy CCD image, with pixel values in keV


%% Input handling

if nargin<1
    error('Need an input argument')
elseif isnumeric(varargin{1})
    img = varargin{1};
else
    error('img input should be numeric array')
end

%defaults
blsig = 0.01874;    %keV. Nov 2007 data gives 0.01874.
fano = 0.14;
epcc = 3.7; %eV here

i=2;
while true
    if i>nargin
        break
    end
    if ~ischar(varargin{i})
        error('Additional arguments must be identified by a string')
    elseif i==nargin    %no following argument
        error('Need argument ID string followed by argument value')
    elseif ~isnumeric(varargin{i+1})
        error('Input argument should be numeric')
    elseif strcmpi(varargin{i},'blsig')
        blsig = varargin{i+1};
    elseif strcmpi(varargin{i},'blfwhm')
        blsig = varargin{i+1} / 2.355;
    elseif strcmpi(varargin{i},'fano')
        fano = varargin{i+1};
    elseif strcmpi(varargin{i},'epcc')
        epcc = varargin{i+1};
    else
        error('Unrecognized argument ID string')
    end
    i = i+2;
end
epcc = epcc / 1e3;  %now in keV

%% Image noise.

u = sqrt(blsig^2 + epcc^2 .* fano .* (img./epcc));    %add stdev's in quadrature to get total sigma of noise
img = img + u.*SampleGauss2(size(img));     %add gaussian noise of magnitude u over image

return

function x = SampleGauss2(s)
%x = SampleGauss2([N,M])
r = rand(s(1),s(2));
x = sqrt(2)*erfinv(r.*2-1);
