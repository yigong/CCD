function Options = HtConstructOptions(nArgs, varArgs)
%
% Parse the varargin from HybridTrack, and assign some builtin values.
%
%   keys        default value
%------------------------------ 
% pixelSizeUm   10.5
% lowThreshold  0.5 [keV]
% dEdxTable     dEdx_ref.mat
% plotflag      false x 15 

Options = struct;

% process input arguments
argNO = 1;
reqArgNO = 1; % required argument number
while argNO <= (nArgs - reqArgNO)
    argInc = 2; % argument increment number
    switch lower(varArgins{argNO})
%        case {'pixelsizeum'}
%            val = varArgs{argNO+1};
%            if isnumeric(val) && isfinite(val) && length(val) == 1
%                Options.pixelSizeUm = val;
%            else
%                error('error in input -- pixelSizeUm');
%            end
%        case {'lowthreshold'}
%            val = varArgs{argNO+1};
%            if isnumeric(val) && isfinite(val) && length(val) == 1
%                Options.lowThreshold = val;
%            else
%                error('error in input -- lowThreshold');
%            end
        case {'plotflags'}
            val = varArgs{argNO+1};
            if islogical(val) || isnumeric(val)
                Options.PlotStyle.multiple = val;
            else
                error('error in input -- plotflags');
            end
        case {'energyt'}
            val = varArgs{argNO+1}
            if isnumeric(val) && isfinite(val) && length(val)==1
                Options.energyT = val;
            else
                error('error in input -- energyT');
        otherwise
            error('error in keyword')
    end
    argNO = argNO + argInc;
end


% set default values -- general
Options.pixelSizeUm = 10.5;
Options.lowThresholdkeV = 0.5; 
Options.dedxTable = HtConstructionDedxTable;
if ~isfield(Options, 'PlotStyle')
    Options.PlotStyle.multiple = false(1,20);
end

% set default values -- ridge following
Options.startDisUm= 30;
Options.ridgeStepLenPix = 0.25; 
Options.cutTotalLenPix = 10;
Options.cutAngleIncDeg = 1;
Options.cutAngleNO = int16(48/Options.cutAngleIncDeg);
Options.piIndices = 180/Options.cutAngleIncDeg;




























