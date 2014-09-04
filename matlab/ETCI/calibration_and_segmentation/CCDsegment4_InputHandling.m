function [imageFull, threshold, opts] = CCDsegment4_InputHandling(nargin, varargin)
%function [imageFull, threshold, PreCalMode, quadrantStraddleVeto, NumPixelsThreshold, NeighborLayers,
%    edgeVeto, subArrayParams, maskFlag, maxSegments, pixelSize, UseBWLabel] = CCDsegment4_InputHandling(nargin, varargin)

%Apparently, varargin gets an extra layer of cell structure when passed as an argument.
% ****NOTE: this is because argument named "varargin" takes a variable number of
%   arguments!!! rename it!
% varargin: 1x1 cell: 1x2 cell
varargin = varargin{1};

%After the first two required arguments (img, th), look for parameter names 
%   and then read the parameter value from the next argument.
%There are some smarter ways to do this, but as it is, each parameter name
%   must be paired with the parameter value.
if nargin < 2
    error('Need at least two input arguments.')
end

%First two arguments are fixed.
imageFull = varargin{1};
threshold = varargin{2};

%initialize options structure.
opts = struct;
i=3;

while i<=nargin     %each argument
    %look at the parameter name, identify it, and make sure that the 
    %   argument value is the correct data type.
    switch lower(varargin{i})
        case 'mode'
            if strcmpi(varargin{i+1},'precal')
                opts.preCalMode = true;
                opts.quadrantStraddleVeto = true;    %veto tracks crossing quadrant edges, if we are performing a calibration
            elseif strcmpi(curarg,'postcal')
                opts.preCalMode = false;
            else
                error('Could not understand given mode. Mode should be ''precal'' or ''postcal'' ')
            end
            i = i+2;
        case 'pixelthreshold'
            if isnumeric(varargin{i+1}) && length(varargin{i+1})==1
                opts.numPixelsThreshold = varargin{i+1};
                i = i + 2;
            else
                error('Could not understand the given pixel threshold.')
            end
        case 'neighborflag'
            if nargin > i && (isnumeric(varargin{i+1}) || islogical(varargin{i+1})) && length(varargin{i+1})==1
                opts.neighborLayers = +varargin{i+1};
                i=i+2;
            else
                opts.neighborLayers = 1;
                i=i+1;
            end
        case 'neighborlayers'
            if nargin > i && isnumeric(varargin{i+1}) && length(varargin{i+1})==1
                opts.neighborLayers = varargin{i+1};
                i=i+2;
            else
                error('Could not understand NeighborLayers input.')
            end
        case 'edgeveto'
            if nargin > i && (isnumeric(varargin{i+1}) || islogical(varargin{i+1})) && length(curarg)==1
                opts.edgeVeto = logical(curarg);
                i=i+2;
            else
                opts.edgeVeto = true;
                i=i+1;
            end
        case 'subarray'
            s = size(varargin{i+1});
            if isnumeric(varargin{i+1}) && all(s==[1,4])
                opts.subArrayParams = varargin{i+1};
            elseif isnumeric(varargin{i+1}) && all(s==[4,1])   %column format?
                opts.subArrayParams = varargin{i+1}';
            elseif isempty(varargin{i+1})  %no subarray
                opts.subArrayParams = [];    %use an empty, instead of a separate flag
            else
                error('Could not understand the subarray parameters.')
            end
            i=i+2;
        case 'maskflag'
            if nargin>i && (isnumeric(varargin{i+1}) || islogical(varargin{i+1})) && length(varargin{i+1})==1
                opts.maskFlag = logical(varargin{i+1});
                i=i+2;
            else
                opts.maskFlag = true;
                i=i+1;
            end
        case 'maxsegments'
            if isnumeric(varargin{i+1}) || islogical(varargin{i+1})
                opts.maxSegments = logical(varargin{i+1});
                i=i+2;
            else
                error('Could not understand the maxSegments argument.')
            end
        case 'pixelsize'
            if isnumeric(varargin{i+1}) && length(varargin{i+1})==1
                opts.pixelSize = varargin{i+1};
                i=i+2;
            else
                error('Could not understand the pixelSize argument.')
            end
        case 'pixsize'  %alternative label to pixelsize
            if isnumeric(varargin{i+1}) && length(varargin{i+1})==1
                opts.pixelSize = varargin{i+1};
                i=i+2;
            else
                error('Could not understand the pixelSize argument.')
            end
        case 'usebwlabel'
            if nargin > i && (isnumeric(varargin{i+1}) || islogical(varargin{i+1})) && length(varargin{i+1})==1
                opts.useBWLabel = logical(varargin{i+1});
                i=i+2;
            else
                opts.useBWLabel = true;
                i=i+1;
            end
        case 'usesmoothing'
            if nargin > i && (isnumeric(varargin{i+1}) || islogical(varargin{i+1})) && length(varargin{i+1})==1
                opts.useSmoothing = logical(varargin{i+1});
                i=i+2;
            else
                opts.useSmoothing = true;
                i=i+1;
            end
        case 'smoothing'
            if nargin > i && (isnumeric(varargin{i+1}) || islogical(varargin{i+1})) && length(varargin{i+1})==1
                opts.useSmoothing = logical(varargin{i+1});
                i=i+2;
            else
                opts.useSmoothing = true;
                i=i+1;
            end
        case 'smoothingkernel'
            if nargin > i && isnumeric(varargin{i+1})
                opts.smoothingKernel = varargin{i+1};
                %assume we want to use this kernel. Unless otherwise specified.
                if ~isfield(opts,'useSmoothing')
                    opts.useSmoothing = true;
                end
                i=i+2;
            else
                error('Could not understand smoothing kernel argument.')
            end
        otherwise
            error('Input argument not recognized.')
    end
end

%Defaults.
if ~isfield(opts,'preCalMode')
    opts.preCalMode = false;
end
if ~isfield(opts,'useBWLabel')
    opts.useBWLabel = false; %use BWConnComp instead
end
if ~isfield(opts,'edgeVeto')
    opts.edgeVeto = true; %also need to look at this again, i think it is written for an old CCDcal code. #futurework
end
if ~isfield(opts,'subArrayParams')
    opts.subArrayParams = [];    %use empty instead of a separate subarray flag
end
if ~isfield(opts,'maskFlag')
    opts.maskFlag = false;
end
if ~isfield(opts,'maxSegments')
    opts.maxSegments = 500;
end
if ~isfield(opts,'pixelSize')
    opts.pixelSize = 10.5;
end
if ~isfield(opts,'NeighborLayers')
	opts.neighborLayers = max(1, round(10.5/opts.pixelSize));
end
if ~isfield(opts,'NumPixelsThreshold')
    opts.numPixelsThreshold = 4; %function of pixel size?
end
if ~isfield(opts,'useSmoothing')
    if opts.pixelSize >= 7.5  %something between 10 and 5
        opts.useSmoothing = false;
    else
        opts.useSmoothing = true;
    end
end
if opts.useSmoothing && ~isfield(opts,'smoothingKernel')
    %get a gaussian-ish kernel with FWHM 10.5 um. 
    %   FWHM here should be somewhat less than the average diffusion, I guess.
    switch opts.pixelSize
        %include special cases to avoid computing time
        case 5
            %fspecial('gaussian', 3, 10.5/5/2.355)
            opts.smoothingKernel = [...
                0.066584897551308   0.124870700849682   0.066584897551308;
                0.124870700849682   0.234177606396040   0.124870700849682;
                0.066584897551308   0.124870700849682   0.066584897551308];
        case 2.5
            %fspecial('gaussian', 5, 10.5/2.5/2.355)
            opts.smoothingKernel = [...
                0.019946887361785   0.031966033134603   0.037407608875818   0.031966033134603   0.019946887361785;
                0.031966033134603   0.051227404849151   0.059947842644442   0.051227404849151   0.031966033134603;
                0.037407608875818   0.059947842644442   0.070152759998388   0.059947842644442   0.037407608875818;
                0.031966033134603   0.051227404849151   0.059947842644442   0.051227404849151   0.031966033134603;
                0.019946887361785   0.031966033134603   0.037407608875818   0.031966033134603   0.019946887361785];
        otherwise
            %hsize is the size of the kernel
            hsize = 10.5/opts.pixelSize;
            %   round up to an odd integer. 
            %   This gives 3x3 for pixelSize >= 3.5um, 5x5 for pixelSize >= 2.1um
            hsize = 1+2*ceil(hsize/2 - 0.5);
            %sigma of gaussian function
            sigma = 10.5/opts.pixelSize/2.355;
            opts.smoothingKernel = fspecial('gaussian',hsize,sigma);
    end
end

opts.quadrantStraddleVeto = opts.preCalMode;
%currently, these are always identical... could be a separate option I suppose.  #futurework

