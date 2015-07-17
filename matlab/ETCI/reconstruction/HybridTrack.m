function Output = HybridTrack(originalImageKev,varargin)
% function Output = HybridTrack(originalImageKev)
% function Output = HybridTrack(originalImageKev, 'key', value, ...)
% 
% Measure the initial electron direction from a CCD image.
% 
% Inputs:
%   originalImageKev: A 2D numeric array of values corresponding to energy
%       deposited in CCD pixels. Units keV.
%     **If pixels are 10.5um and no special operation or plotting is desired,
%       no key/value pairs are needed.
%   key/value: The following key/value pairs are supported, in the format
%       e.g. HybridTrack(img, 'pixelSizeUm', 5, 'lowThresholdKev', 0.3);
%       Defaults given in []. Key names are case-insensitive.
%     pixelSizeUm: Associate a different pixel size with the input image.
%       [10.5]
%     lowThresholdKev: Use a different threshold value for the binary image 
%       step. [0.5 for 10.5um pixels, and scales with pixel area]
%     dedxTable: Use a different reference table of dE/dx vs. energy, for
%       interpolating values to use in calculating beta.
%       Format: column 1 is energy (keV), column 2 is dEdx (keV/um)
%       [uses some rough values grabbed from geant in 2008, a.k.a. dEdx_ref.mat]
%     cheat: Structure of information from geant4, created in the codes
%       Geant4TrackHandling and DiffuseTrack. Not used, but passed through to
%       output structure for calculating algorithm errors afterwards.
%     oneplot: A single logical value. (You can just put this key as a flag,
%       without having to supply a value, and it will be interpreted as 'true'.)
%       Indicates that a plot should be generated showing the ridge following
%       points, best cuts, and step directions.
%     plotflags: Generate plots illustrating any of various steps in the
%       algorithm process. 1x15 logical array, with any true's generating a plot
%       associated with that index. (Up to 15 plots.) I'm not listing them all
%       here at this time. Default all false.
%     insetPix: View coordinates for some of the plots in plotflags. 
%       assign by: insetPix = [xlim, ylim];
%     nameString: For saving plots to file. File names end up in the format:
%       ['Track_', nameString, '_plotSpecificSuffix.eps'] and *.fig and *.png.
% 
% Outputs:
%   Output: structure containing a whole bunch of stuff. The most user-friendly
%       and useful fields are as follows:
%     alpha: the estimated direction in the plane of the pixels. Units degrees.
%       *alpha = 0 is in the positive x direction, where x is the first array
%       dimension. In an image array in Matlab variable editor, this is down; in
%       a surface plot or pseudocolor plot in Matlab, this is up.
%       *alpha = 90 is in the positive y direction, where y is the second array
%       dimension. In an image array in Matlab variable editor, this is to the
%       right; in a surface plot or pseudocolor plot in Matlab, this is also to
%       the right.
%     beta: the estimated direction, as the magnitude of the out-of-plane angle.
%       Units degrees. This is a rather poor estimate right now.
%     img: the input image with a one-pixel buffer around the outside. This is
%       the image to which any reference coordinates are given in.
%     Etot: Etot = sum(originalImageKev(:))
%     

% follows version 1c



%The lowthresh should cut out some diffusion,
% so that the thinned track is closer to the actual track,
% and the correct end can be identified more accurately.
%If the lowthresh is too high, the initial end may be cut off.

%%%%%%%%%%%%%%%
%%%% input %%%%
%%%%%%%%%%%%%%%

Options = HtConstructOptions(nargin,varargin);
% this includes the ridge-following options, in a subroutine.

% add buffer of zeros around image
[trackEnergy, preparedImageKev, Options] = ...
    HtPrepareImage(originalImageKev, Options);

%%%%%%%%%%%%%%%
%%%% ends %%%%%
%%%%%%%%%%%%%%%

% low threshold, thinning, identify ends.
EdgeSegments = HtChooseInitialEnd(preparedImageKev, Options);

% exception
if isnan(EdgeSegments.chosenIndex)
    % no end found
    % exit unsuccessfully
    Output.img = preparedImageKev;
    Output.Etot = sum(preparedImageKev(:));
    Output.ends = 0;
    Output.err = 'No ends found';
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ridge following %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% ridge following
Ridge = HtRidgeFollow(EdgeSegments, Options, preparedImageKev);

% exception
if isfield(Ridge,'err') && strcmpi(Ridge.err,'infinite loop')
    Output.img = preparedImageKev;
    Output.Etot = sum(preparedImageKev(:));
    Output.err = 'Infinite loop';
    Output.EdgeSegments = EdgeSegments;
    Output.Ridge = Ridge;
    Output.Options = Options;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% compute direction %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Measurement, Ridge] = HtComputeDirection(trackEnergy, Ridge, Options);

%%%%%%%%%%%%%%%%
%%%% output %%%%
%%%%%%%%%%%%%%%%
Output = ...
    HtSetOutput(preparedImageKev, Options, EdgeSegments, Ridge, Measurement);

%%%%%%%%%%%%%%
%%%% plot %%%%
%%%%%%%%%%%%%%
% plots #13-15
if Options.PlotStyle.multiple(13)
    % all ridge points, measurement highlighted
    [~, Options.PlotStyle.axesPosition] = ...
        HtPlotImage(preparedImageKev, Options.PlotStyle.cmap, ...
        Options.PlotStyle.axesPosition, ...
        Options.PlotStyle.insetPix);
    hold on
    % points
    plot3(Ridge.positionPix(:,2)+0.5, Ridge.positionPix(:,1)+0.5, ...
        Options.PlotStyle.z0*ones(size(Ridge.positionPix(:,1))), ...
        'd','markerfacecolor',Options.PlotStyle.ptColor, ...
        'markeredgecolor',Options.PlotStyle.ptColor, ...
        'markersize',4);
    % highlighted measurement
    measInd = Measurement.indices;
    plot3(Ridge.positionPix(measInd,2)+0.5, Ridge.positionPix(measInd,1)+0.5, ...
        Options.PlotStyle.z0*ones(size(Ridge.positionPix(measInd,1))), ...
        'd','markerfacecolor',Options.PlotStyle.measPtColor, ...
        'markeredgecolor',Options.PlotStyle.measPtColor, ...
        'markersize',4);
    HtSavePlot(gcf, Options, '13_meas_highlighted')
end
if Options.PlotStyle.multiple(14)
    % dE/dx plot, measurement highlighted, dE/dx values marked horizontally
    HtPlotImage(ones(2),'gray',Options.PlotStyle.axesPosition, []);
    hold off
    % measured dE points
    xtmp = cumsum(Ridge.stepLengthPix);
    ytmp = Ridge.dedxKevUm;
    plot(xtmp, ytmp, '*', 'color', Options.PlotStyle.ptColor, ...
        'linewidth',2, 'markersize',8);
    hold on
    % dE/dx horizontal lines
    dedxLineWidth = 2;
    hLegend(1) = plot(xtmp, Measurement.dedxMeasured*ones(size(xtmp)), '-', ...
        'color',Options.PlotStyle.measPtColor, 'linewidth',dedxLineWidth);
    hLegend(2) = plot(xtmp, Measurement.dedxReference*ones(size(xtmp)), '-', ...
        'color','k', 'linewidth',dedxLineWidth);
    if ~isempty(Options.cheat)
        ytmp = Measurement.dedxReference / cosd(Options.cheat.beta);
        hLegend(3) = plot(xtmp, ytmp*ones(size(xtmp)), '-', ...
            'color',Options.PlotStyle.trueArrColor, 'linewidth',dedxLineWidth);
    end
    xlabel('Distance along ridge [pixel lengths]');
    ylabel('dE/dx [keV / {\mu}m]')
    xlim([0, xtmp(end)+Options.positionStepSizePix]);
    hLegend(4) = legend(hLegend, {'Measured \Delta{E}/\Delta{R}', ...
        '(dE/dS)_{calculated}', ...
        '(dE/dS)_{calculated} / cos(\beta)'}, ...
        'Location','SouthEast');
    set(hLegend(4),'fontsize',16);
    FixAxesMargins(gca);
    HtSavePlot(gcf, Options, '14_dEdx_meas');
end
if Options.PlotStyle.multiple(15)
    % ridge points, measurement highlighted, with arrow(s)
    [~, Options.PlotStyle.axesPosition] = ...
        HtPlotImage(preparedImageKev, Options.PlotStyle.cmap, ...
        Options.PlotStyle.axesPosition, ...
        Options.PlotStyle.insetPix);
    hold on
    % points
    plot3(Ridge.positionPix(:,2)+0.5, Ridge.positionPix(:,1)+0.5, ...
        Options.PlotStyle.z0*ones(size(Ridge.positionPix(:,1))), ...
        'd','markerfacecolor',Options.PlotStyle.ptColor, ...
        'markeredgecolor',Options.PlotStyle.ptColor, ...
        'markersize',4);
    % measurement highlighted
    measInd = Measurement.indices;
    plot3(Ridge.positionPix(measInd,2)+0.5, Ridge.positionPix(measInd,1)+0.5, ...
        Options.PlotStyle.z0*ones(size(Ridge.positionPix(measInd,1))), ...
        'd','markerfacecolor',Options.PlotStyle.measPtColor, ...
        'markeredgecolor',Options.PlotStyle.measPtColor, ...
        'markersize',4);
    % measured arrow
    arr_length1 = 4;    %pixels. main stem.
    arr_length2 = .75;  %pixels. Side prongs.
    arrowLineWidth = 2;
    xTmp0 = Ridge.positionPix(measInd(1),1)+0.5;
    yTmp0 = Ridge.positionPix(measInd(1),2)+0.5;
    %   draw arrow from base to tip to two prongs
    a = Measurement.alphaDegrees;
    xMeasArrow = xTmp0 + [0, ...
        cosd(a)*arr_length1, ...
        cosd(a)*arr_length1 - cosd(a+45)*arr_length2, ...
        cosd(a)*arr_length1, ...
        cosd(a)*arr_length1 - cosd(a-45)*arr_length2];
    yMeasArrow = yTmp0 + [0, ...
        sind(a)*arr_length1, ...
        sind(a)*arr_length1 - sind(a+45)*arr_length2, ...
        sind(a)*arr_length1, ...
        sind(a)*arr_length1 - sind(a-45)*arr_length2];
    plot3(yMeasArrow, xMeasArrow, ...
        Options.PlotStyle.z0*ones(size(yMeasArrow)), ...
        '-', 'color',Options.PlotStyle.measArrColor, ...
        'linewidth',arrowLineWidth);
    if ~isempty(Options.cheat)
        % true arrow
        a = Options.cheat.alpha;
        xTrueArrow = xTmp0 + [0, ...
            cosd(a)*arr_length1, ...
            cosd(a)*arr_length1 - cosd(a+45)*arr_length2, ...
            cosd(a)*arr_length1, ...
            cosd(a)*arr_length1 - cosd(a-45)*arr_length2];
        yTrueArrow = yTmp0 + [0, ...
            sind(a)*arr_length1, ...
            sind(a)*arr_length1 - sind(a+45)*arr_length2, ...
            sind(a)*arr_length1, ...
            sind(a)*arr_length1 - sind(a-45)*arr_length2];
        plot3(yTrueArrow, xTrueArrow, ...
            Options.PlotStyle.z0*ones(size(yTrueArrow)), ...
            '-', 'color',Options.PlotStyle.trueArrColor, ...
            'linewidth',arrowLineWidth);
    end
    HtSavePlot(gcf, Options, '15_alpha_arrows');
end



function Options = HtConstructOptions(nArgs, varargs)
% function InputOptions = HtConstructOptions(nArgs, varargs)
%
% Parse the varargin from HybridTrack, and assign some builtin values.
% 
% key-value pairs
%
% keys (alternate strings) [default value]:
%   pixelSizeUm (pixelsize, pixsize, pixelpitch, pixel, pitch)
%       The pixel size in the image, in microns.. Used to scale many other 
%       parameters.
%       [10.5]
%   lowThreshold (lt, lowthresh, threshold)
%       The threshold to apply to make the binary image. 
%       [0.5 keV * (pixelSizeUm/10.5)^2]
%   dEdxTable (dEdx_table, dEdx_ref, dEref, dE_ref, dEdx)
%       The table of linear stopping powers, as a function of energy, to use for
%       calculating beta. First column, energy in keV, second column, 
%       stopping power in keV/um.
%       [defined from 2009(?) values in dEdx_ref.mat]
%   plotflag
%       ....define how to approach this...
%       [false]
%   cheat
%       Structure of information from Geant4TrackHandling, providing true values
%       for electron track.
%       []
% 
% Also, call HtDefineDefaultRidgeFollowingOptions

Options = struct;

%first, parse input arguments
argumentNo = 1;
requiredArgs = 1;   %because varargin includes required arguments too
while argumentNo <= (nArgs - requiredArgs);
    argumentIncrement = 2;  %by default, increment the argumentNo by 2 each time
    switch lower(varargs{argumentNo})
        case {'pixelsizeum','pixelsize','pixsize','pixelpitch','pixel','pitch'}
            thisValue = varargs{argumentNo+1};
            if isnumeric(thisValue) && isfinite(thisValue) && ...
                    length(thisValue)==1
                Options.pixelSizeUm = thisValue;
            else
                error('pixelSizeUm must be a single value, numeric and finite');
            end
        case {'lowthresh','lt','lowthreshold','threshold','lowthresholdkev'}
            thisValue = varargs{argumentNo+1};
            if isnumeric(thisValue) && isfinite(thisValue) && ...
                    length(thisValue)==1
                Options.lowThreshold = thisValue;
            else
                error('lowThresholdKev must be a single value, numeric and finite')
            end
        case {'dedx_ref','dedx','deref','de_ref','dedxtable','dedx_table'}
            thisValue = varargs{argumentNo+1};
            if isnumeric(thisValue) && all(isfinite(thisValue(:))) && ...
                    size(thisValue,2)>=2    %two columns of values
                Options.dedxTable = thisValue;
            else
                error('dedxTable should be two columns, numeric and finite')
            end
        case {'plotflag', 'plotflags'}
            if argumentNo+1 > (nArgs - requiredArgs) || ...
                    ischar(varargs{argumentNo+1})
                % no value is supplied. this means, set PlotStyle.single = true
                Options.PlotStyle.single = true;
                Options.PlotStyle.multiple = false(1,20);
                argumentIncrement = 1;
            elseif islogical(varargs{argumentNo+1}) && ...
                    length(varargs{argumentNo+1})==1
                Options.PlotStyle.single = varargs{argumentNo+1};
                Options.PlotStyle.multiple = false(1,20);
            elseif islogical(varargs{argumentNo+1}) && ...
                    length(varargs{argumentNo+1}) > 1
                Options.PlotStyle.single = false;
                Options.PlotStyle.multiple = varargs{argumentNo+1};
            elseif isnumeric(varargs{argumentNo+1})
                Options.PlotStyle.single = false;
                Options.PlotStyle.multiple = false(1,20);
                Options.PlotStyle.multiple(varargs{argumentNo+1}) = true;
            else
                error('plotflag should be a boolean or numeric value')
            end
        case {'insetpix','inset','xylim'}
            thisArg = varargs{argumentNo+1};
            if isnumeric(thisArg) && ...
                    length(thisArg)==4
                Options.PlotStyle.insetPix = thisArg(1:4);
            elseif ~isempty(thisArg)
                error('insetpix should be a 1x4 numeric array')
            end
        case {'namestring','savename','name','string','plotname'}
            thisArg = varargs{argumentNo+1};
            if ischar(thisArg) && isempty(strfind(thisArg,'*'))
                Options.PlotStyle.nameString = thisArg;
                % Track_*_descriptor.{fig|png|eps}
            elseif ischar(thisArg)
                error('Don''t put special characters in the filename! (*)')
            end
        case {'cheat'}
            thisValue = varargs{argumentNo+1};
            if isstruct(thisValue) % && length(thisValue)==1
                Options.cheat = thisValue;
            else
                error('cheat should be a structure')
            end
        case {'oneplot', 'plot'}
            if argumentNo+1 > (nArgs - requiredArgs) || ischar(varargs{argumentNo+1})
                Options.PlotStyle.single = true;
                Options.PlotStyle.multiple = false(1,20);
                argumentIncrement = 1;
            elseif islogical(varargs{argumentNo+1}) || ...
                    isnumeric(varargs{argumentNo+1})
                Options.PlotStyle.single = varargs{argumentNo+1};
                Options.PlotStyle.multiple = false(1,20);
            else
                error('unrecognized plot style')
            end
        case {'einit', 'energyt'}
            thisValue = varargs{argumentNo+1};
            if isnumeric(thisValue) && isfinite(thisValue) ...
                    && length(thisValue)==1
                Options.energyT = thisValue; % in keV
            else
                error('unrecognized energy')
            end
        otherwise
            error('Unrecognized input argument')
    end
    argumentNo = argumentNo + argumentIncrement;
end

%now, assign default values where needed
if ~isfield(Options,'pixelSizeUm')
    Options.pixelSizeUm = 10.5;
end
if ~isfield(Options,'lowThresholdKev')
    Options.lowThresholdKev = 0.5 * (Options.pixelSizeUm / 10.5)^2;
end
if ~isfield(Options,'dedxTable')
    Options.dedxTable = HtConstructDedxTable;
end
if ~isfield(Options,'plotflag')
    Options.plotflag = false;
end
if ~isfield(Options,'cheat')
    Options.cheat = [];
end
if ~isfield(Options,'PlotStyle')
    Options.PlotStyle.single = false;
    Options.PlotStyle.multiple = false(1,20);
end

%assign additional options which are not actually inputs (yet)
Options = HtDefineDefaultRidgeFollowingOptions(Options);
Options = HtDefineMeasurementOptions(Options);

Options.ridgeStartingDistanceFromTrackEndUm = 40;

Options = HtDefinePlotOptions(Options);



function dedxTable = HtConstructDedxTable
% function dedxTable = constructDedxTable
% 
% Create a lookup table for dEdx_reference.
% 
% Column 1: electron energy [keV]
% Column 2: dEdx_reference [keV/um]
% 
% Values from dEdx_ref.mat, from 2009(?) work.

% TODO: Revisit this, preferably using theory.

dedxTable(:,1) = 50:25:1400;    %keV
dedxTable(:,2) = [ 1.2114       %keV / um
    0.9185
    0.7200
    0.5800
    0.5300
    0.4400
    0.3700
    0.3600
    0.3200
    0.3100
    0.2900
    0.2900
    0.2700
    0.2700
    0.2600
    0.2500
    0.2500
    0.2400
    0.2400
    0.2400
    0.2300
    0.2300
    0.2300
    0.2300
    0.2200
    0.2200
    0.2200
    0.2200
    0.2200
    0.2200
    0.2200
    0.2200
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100
    0.2100];



function Options = HtDefineDefaultRidgeFollowingOptions(Options)
% function Options = HtDefineDefaultRidgeFollowingOptions(Options)
% 
% Assign fixed parameters used for the ridge-following.
% 
% Input:
%   Options structure defined in HtConstructOptions
%   must include field pixelSizeUm
% 
% Output structure includes the following additional fields:
%   positionStepSizePix
%   cutSamplingIntervalPixPix
%   cutTotalLengthPix
%   cutAngleIncrementDegrees
%   cutAngleIncrementRadians
%   searchAngleIndices
%   piIndices
%   cutAngleDegrees
%   cutAngleRadians
%   cutInterpolationMethod
%   cutDistanceFromCenterPix
%   cutCoordinatesPix
%   trackEndLowThresholdKev

Options.positionStepSizePix = 0.25;
Options.cutSamplingIntervalPix = 0.25;

% (smaller cutTotalLengthPix will make code run faster)
Options.cutTotalLengthPix = 105 / Options.pixelSizeUm;

% This needs to be a factor of 45
% Smaller cutAngleIncrementDegrees might give a more accurate measurement of 
% alpha.
% (larger cutAngleIncrementDegrees will make code run faster)
Options.cutAngleIncrementDegrees = 3;

% cutAngleIncrementDegrees is used as the base unit for all the other anglular
% variables. So variables ending in "Indices" are in this unit.

% This needs to be a multiple of 2. The units are, the number of indices of
%  cutAngleIncrement.
% The maximum angular change in a single step is searchAngleIndices/2.
% TODO: this could be brought down for smaller pixsize, depending on PSS
% Larger searchAngleIndices will allow the ridge-following to make a tighter 
%  turn. However at 4 steps per pixel length, it does not need to make a very 
%  tight turn. Regardless, searchAngleIndices/2 should always be much less than 
%  90 degrees or else there is more potential for losing the ridge and walking 
%  off an elbow or something.
% (smaller searchAngleIndices will make code run faster)
Options.searchAngleIndices = 48/Options.cutAngleIncrementDegrees;  

Options.twoPiIndices = 360/Options.cutAngleIncrementDegrees;
Options.piIndices = 180/Options.cutAngleIncrementDegrees;

Options.cutAngleIncrementRadians = Options.cutAngleIncrementDegrees * pi/180;

Options.cutAngleRadians(1:Options.twoPiIndices) = ...
    Options.cutAngleIncrementRadians : ...
    Options.cutAngleIncrementRadians : ...
    Options.twoPiIndices*Options.cutAngleIncrementRadians;
Options.cutAngleDegrees(1:Options.twoPiIndices) = ...
    Options.cutAngleIncrementDegrees : ...
    Options.cutAngleIncrementDegrees : ...
    Options.twoPiIndices*Options.cutAngleIncrementDegrees;

% TODO: Options.cutLowThresholdKev, for truncating cut width

Options.cutInterpolationMethod = 'linear';
% other options include 'cubic' and 'spline', which would be slower.

% Define x,y of an angled cut for a step direction of 0 degrees, then rotate to
% all possible step angles.
cut0y = -Options.cutTotalLengthPix/2 : ...
    Options.cutSamplingIntervalPix : ...
    Options.cutTotalLengthPix/2;
cut0x = zeros(size(cut0y));

% distance from center to a point. This is used for measuring width, 
%   width = sum_i (d_i * E_i)
Options.cutDistanceFromCenterPix = abs(cut0y);
Options.cutDistanceCoordinatePix = cut0y;

Options.cutCoordinatesPix = cell(1,Options.twoPiIndices);
for angleIndex=1:Options.twoPiIndices
    % rotation matrix
    R = [cos(Options.cutAngleRadians(angleIndex)), ...
         sin(Options.cutAngleRadians(angleIndex)); ...
        -sin(Options.cutAngleRadians(angleIndex)), ...
         cos(Options.cutAngleRadians(angleIndex))];
    Options.cutCoordinatesPix{angleIndex} = [cut0x',cut0y']*R;
end

baseCutLowThresholdKev = 0.05;
Options.cutLowThresholdKev = baseCutLowThresholdKev * ...
    (Options.pixelSizeUm / 10.5).^2;
% cut points below this threshold, and beyond, are ignored

baseTrackEndLowThresholdKev = 0.1;
Options.trackEndLowThresholdKev = baseTrackEndLowThresholdKev * ...
    (Options.pixelSizeUm / 10.5)^2;     % scale with pixel area.

Options.infiniteLoopThresholdPix = Options.positionStepSizePix / 2;



function Options = HtDefineMeasurementOptions(Options)
% function Options = HtDefineMeasurementOptions(Options)

minimumWidthMeasurementLengthUm = 30;   % roughly defined by diffusion
preferredWidthMeasurementLengthPix = 2;
Options.widthMeasurementLengthPix = ...
    max(preferredWidthMeasurementLengthPix, ...
        minimumWidthMeasurementLengthUm / Options.pixelSizeUm);

Options.initialBetaGuessDegrees = 0;
Options.shouldShortenMeasurementLength = true;
Options.measurementFunctionHandle = @median;



function Options = HtDefinePlotOptions(Options)
% function Options = HtDefinePlotOptions(Options)

Options.PlotStyle.axesPosition = [];    % gets set in first plot call
if any(Options.PlotStyle.multiple)
    load cmaps.mat;
    Options.PlotStyle.cmap = cmaphotlog;
    if ~isfield(Options.PlotStyle,'insetPix')
        Options.PlotStyle.insetPix = [];
    end
    if ~isfield(Options.PlotStyle,'nameString')
        Options.PlotStyle.nameString = [];
    end
end

Options.PlotStyle.z0 = 50;

% colors
Options.PlotStyle.cutColor = 'y';
Options.PlotStyle.bestCutColor = 'g';
Options.PlotStyle.otherCutColor = 'c';
Options.PlotStyle.ptColor = [0,0.6,1];
Options.PlotStyle.measPtColor = 'g';
Options.PlotStyle.measArrColor = 'g';
Options.PlotStyle.trueArrColor = [0.5,0,1];



function [trackEnergyKev, newImageKev, Options] = ...
    HtPrepareImage(originalImageKev, Options)
% function [trackEnergy, newImage, Options] = ...
%     HtPrepareImage(originalImageKev, Options)
% 
% Add buffer around track image.
% 
% Inputs:
%   Options structure defined in HtConstructOptions
%   must include field cutTotalLengthPix

% imageEdgeBuffer does not need to handle the cutTotalLength/2 in any direction.
%   it only needs to handle the ridge points going off the edge.

% TODO: but here, we operate under a misunderstanding. 
imageEdgeBufferUm = 0.55 * Options.cutTotalLengthPix;
imageEdgeBufferPix = ceil(imageEdgeBufferUm / Options.pixelSizeUm);

newImageSize = [size(originalImageKev,1) + 2*imageEdgeBufferPix, ...
                size(originalImageKev,2) + 2*imageEdgeBufferPix];
newImageKev = zeros(newImageSize);

xIndicesOfOriginalImage = imageEdgeBufferPix + (1:size(originalImageKev,1));
yIndicesOfOriginalImage = imageEdgeBufferPix + (1:size(originalImageKev,2));
newImageKev(xIndicesOfOriginalImage, yIndicesOfOriginalImage) = originalImageKev;
if isfield(Options, 'energyT')
    trackEnergyKev = Options.energyT;
    Options.energyM = sum(newImageKev(:));
else
    trackEnergyKev = sum(newImageKev(:));
    
end



function EdgeSegments = HtChooseInitialEnd(imgKev, Options)
% function EdgeSegments = HtChooseInitialEnd(imgKev, Options)
% 
% Analyze CCD image to decide what end to measure.
% 
% Inputs:
%   imgKev: CCD image.
%   Options: structure from HtConstructOptions.
% 
% Outputs: EdgeSegments structure with the following fields:
%   energiesKev: list of energies measured at each end
%   coordinatesPix: Nx2 list of x,y coordinates for each end pixel
%   chosenIndex: index for above two fields, of the selected end
%   startRidgePix: x,y coordinates to start the ridge following
%   startDirectionIndices: direction to start ridge following, based on thinned
%       track image, in units of cutAngleIncrementDegrees

% first: locate all ends
endLinearIndices = [];
connectivity = ones(3);     %8-connectivity

% if we are unsuccessful in finding ends, increase threshold to break loops.
% Normally this loop will only run once.
lowThresholdKev = Options.lowThresholdKev;    %initial value only
while isempty(endLinearIndices) && lowThresholdKev <= 10*Options.lowThresholdKev;
    binaryImage = +(imgKev > lowThresholdKev);
    thinnedImage = +bwmorph(binaryImage,'thin',inf);
    nNeighborsImage = conv2(thinnedImage,connectivity) - 1;
    nNeighborsImage = nNeighborsImage(2:end-1,2:end-1); %reduce to original size
    nNeighborsImage = nNeighborsImage .* thinnedImage;  %only pixels in the track
    
    endImage = (nNeighborsImage==1);
    endLinearIndices = find(endImage);
    
    %if we have not found any ends, increment by the original threshold.
    if isempty(endLinearIndices)
        lowThresholdKev = lowThresholdKev + Options.lowThresholdKev;
    end
end

if ~isempty(endLinearIndices)
    lowThresholdUsed = lowThresholdKev;
else
    %error: no ends found
    EdgeSegments.energiesKev = [];
    EdgeSegments.coordinatesPix = zeros(0,2);
    EdgeSegments.chosenIndex = nan;
    EdgeSegments.startRidgePix = zeros(0,2);
    EdgeSegments.startDirectionIndices = nan;
    EdgeSegments.lowThresholdUsed = nan;
    return
end

% measure energies
[psf,psfArrayRadiusPix] = HtConstructPsf(Options.pixelSizeUm);
energySumImage = conv2(imgKev,psf);
energySumImage = energySumImage(1+psfArrayRadiusPix:end-psfArrayRadiusPix, ...
                                1+psfArrayRadiusPix:end-psfArrayRadiusPix);
EdgeSegments.energiesKev = energySumImage(endLinearIndices);
[EdgeSegments.coordinatesPix(:,1), EdgeSegments.coordinatesPix(:,2)] ...
    = ind2sub(size(imgKev),endLinearIndices);
[~,EdgeSegments.chosenIndex] = min(EdgeSegments.energiesKev);
EdgeSegments.chosenEnd = EdgeSegments.coordinatesPix(EdgeSegments.chosenIndex, :);
% Walk back up the track to get starting location. For each step, I count
% neighbors to see if we're at an intersection. If not, I use a find operation
% to get the position of the next pixel
nStepsPixels = ceil(Options.ridgeStartingDistanceFromTrackEndUm / ...
                    Options.pixelSizeUm);
thisEnd = EdgeSegments.coordinatesPix(EdgeSegments.chosenIndex, :);
imageTemp = thinnedImage;
thisXY = thisEnd;
% keep track of the direction of every step.
% step #1 is from pixel #1 (end pixel) to pixel #2.
xStep = nan(1,nStepsPixels);    
yStep = nan(1,nStepsPixels);
for stepNo = 1:nStepsPixels
    imageTemp(thisXY(1), thisXY(2)) = 0;
    neighbors = imageTemp(thisXY(1)-1:thisXY(1)+1, thisXY(2)-1:thisXY(2)+1);
    nNeighbors = sum(neighbors(:));
    if nNeighbors == 1
        % take this step
        % by definition, this must happen for the first step, since we start at
        % an end.
        [xStep(stepNo), yStep(stepNo)] = find(neighbors);
        xStep(stepNo) = xStep(stepNo) - 2;  %from [1 2 3] to [-1 0 1]
        yStep(stepNo) = yStep(stepNo) - 2;  %from [1 2 3] to [-1 0 1]
        thisXY = thisXY + [xStep(stepNo), yStep(stepNo)];
        lastStep = [xStep(stepNo), yStep(stepNo)];  %for direction
    elseif nNeighbors > 1
        % at an intersection.
        % take a step back.
        thisXY = thisXY - lastStep;
        % lastStep can remain as is, pointing away from the intersection.
        break
    elseif nNeighbors == 0
        % end of track. keep thisXY
        break
    end
end
EdgeSegments.startRidgePix = thisXY;
startDirectionDegrees = 180/pi * atan2(-lastStep(2), -lastStep(1));
if startDirectionDegrees > 0    %atan2 returns a value on [-pi, pi]
    EdgeSegments.startDirectionIndices = round(startDirectionDegrees / ...
        Options.cutAngleIncrementDegrees);
else
    EdgeSegments.startDirectionIndices = round((startDirectionDegrees + 360) ...
        / Options.cutAngleIncrementDegrees);
end

EdgeSegments.lowThresholdUsed = lowThresholdUsed;
EdgeSegments.thinnedImage = thinnedImage;

if Options.PlotStyle.multiple(1)
    % original image
    [~, Options.PlotStyle.axesPosition] = ...
        HtPlotImage(imgKev, Options.PlotStyle.cmap, ...
        Options.PlotStyle.axesPosition, []);
    HtSavePlot(gcf, Options, '01_original');
end
if Options.PlotStyle.multiple(2)
    % standard binary image
    binImg = imgKev > Options.lowThresholdKev;
    [~, Options.PlotStyle.axesPosition] = ...
        HtPlotImage(binImg, 'gray', ...
        Options.PlotStyle.axesPosition, []);
    HtSavePlot(gcf, Options, '02_binary');
end
if Options.PlotStyle.multiple(3)
    % high-threshold binary image
    binImg = imgKev > 3*Options.lowThresholdKev;
    [~, Options.PlotStyle.axesPosition] = ...
        HtPlotImage(binImg, 'gray', ...
        Options.PlotStyle.axesPosition, []);
    HtSavePlot(gcf, Options, '03_binary_high');
end
if Options.PlotStyle.multiple(4)
    % low-threshold binary image
    binImg = imgKev > 0.1*Options.lowThresholdKev;
    [~, Options.PlotStyle.axesPosition] = ...
        HtPlotImage(binImg, 'gray', ...
        Options.PlotStyle.axesPosition, []);
    HtSavePlot(gcf, Options, '04_binary_low');
end
if Options.PlotStyle.multiple(5)
    % thinned image
    binImg = EdgeSegments.thinnedImage;
    [~, Options.PlotStyle.axesPosition] = ...
        HtPlotImage(binImg, 'gray', ...
        Options.PlotStyle.axesPosition, []);
    HtSavePlot(gcf, Options, '05_thinned');
end
if Options.PlotStyle.multiple(6)
    % thinned image, ends highlighted
    binImg = +EdgeSegments.thinnedImage;
    binImg(endLinearIndices) = 2;
    [~, Options.PlotStyle.axesPosition] = ...
        HtPlotImage(binImg, 'gray', ...
        Options.PlotStyle.axesPosition, []);
    HtSavePlot(gcf, Options, '06_thinned_highlighted');
end
if Options.PlotStyle.multiple(7)
    % image with ends circled
    [h, Options.PlotStyle.axesPosition] = ...
        HtPlotImage(imgKev, Options.PlotStyle.cmap, ...
        Options.PlotStyle.axesPosition, []);
    for i=1:size(EdgeSegments.coordinatesPix,1)
        x0tmp = EdgeSegments.coordinatesPix(i,1);
        y0tmp = EdgeSegments.coordinatesPix(i,2);
        x0tmp = x0tmp + 0.5;
        y0tmp = y0tmp + 0.5;
        disp(['End #',num2str(i),' at (',num2str(x0tmp),', ',num2str(y0tmp),...
            '): energy = ',num2str(EdgeSegments.energiesKev(i))])
        %draw the outline, starting at lower left inside corner
        xtmp = [x0tmp-1.5, x0tmp-2.5, x0tmp-2.5, x0tmp-1.5, x0tmp-1.5, x0tmp+1.5, x0tmp+1.5, ...
            x0tmp+2.5, x0tmp+2.5, x0tmp+1.5, x0tmp+1.5, x0tmp-1.5, x0tmp-1.5]';
        ytmp = [y0tmp-1.5, y0tmp-1.5, y0tmp+1.5, y0tmp+1.5, y0tmp+2.5, y0tmp+2.5, y0tmp+1.5, ...
            y0tmp+1.5, y0tmp-1.5, y0tmp-1.5, y0tmp-2.5, y0tmp-2.5, y0tmp-1.5]';
        hold on;
        plot3(h, ytmp, xtmp, Options.PlotStyle.z0*ones(size(xtmp)), ...
            'g','linewidth',1);
    end
    HtSavePlot(gcf, Options, '07_ends_circled');
end



function [psf, psfArrayRadiusPix] = HtConstructPsf(pixelSizeUm)
% function [psf, psfArrayRadiusPix] = HtConstructPsf(pixelSizeUm)
% 
% Construct a convolution kernel covering a 25 um radius.
% 
% Inputs:
%   pixelSizeUm: pixel size in microns
% 
% Outputs:
%   psf: the convolution kernel, consisting of 0s and 1s.
%   psfArrayRadiusPix: (size(psf,1)-1)/2

%1.0
%Sum energy from pixels within 25 um of end pixel.
%   These are by no means tested or optimized . . this is something to work on.
%
% Viewing code for circle: (command window)
% pixsize = 10.5; th = 0:360; x = 0.5+25/pixsize*cosd(th); y = 0.5+25/pixsize*sind(th); plot(x,y)
% set(gca,'xtick',-ceil(1+25/pixsize):ceil(1+25/pixsize)); set(gca,'ytick',-ceil(1+25/pixsize):ceil(1+25/pixsize)); 
% 
% This code replicates the hard-coded psf's for 5, 10.5, 20, and 30+ um. 

psfRadiusUm = 25;
psfExactRadiusPix = psfRadiusUm / pixelSizeUm;
psfArrayRadiusPix = ceil(psfExactRadiusPix) - 1;    %...
psfDiameterPix = psfArrayRadiusPix * 2 + 1;  %odd integer
x = -psfArrayRadiusPix : psfArrayRadiusPix;
y = -psfArrayRadiusPix : psfArrayRadiusPix;
psf = zeros(psfDiameterPix);
for xInd = 1:length(x)
    for yInd = 1:length(y)
        if x(xInd)^2 + y(yInd)^2 < psfExactRadiusPix^2
            psf(xInd,yInd) = 1;
        end
    end
end



function [h,axSize] = HtPlotImage(imageToPlot,cmap,axSize,insetPix)
% function [h,axSize] = HtPlotImage(imageToPlot,cmap,axSize,insetPix)
% 
% imageToPlot: numeric or boolean
% cmap: color map (use cmaphotlog from cmaps.mat)
% 
% h: axes handle

axesFontSize = 18;
xSizeScale = 1;
ySizeScale = 1;
figPosition = [851, 388, xSizeScale*560, ySizeScale*420];

grayColorMap = ischar(cmap) && strcmpi(cmap,'gray');
if isnumeric(imageToPlot) && ~grayColorMap
    % normal track image
    h = SurfElectronTrack(imageToPlot,'cmap',cmap);
elseif islogical(imageToPlot) || grayColorMap
    % binary image
    h = SurfElectronTrack(+imageToPlot,'cmap',cmap);
    colorbar('off');
else
    error('imageToPlot should be numeric or logical')
end
if ~isempty(insetPix)
    xlim(h,insetPix(1:2));
    ylim(h,insetPix(3:4));
end
set(h,'fontsize',axesFontSize);
set(get(h,'parent'),'position',figPosition);
drawnow;
if isempty(axSize)
    % first time only
    axSize = get(h,'position');
end
set(h,'position',axSize);
xlabel([]);
ylabel([]);
set(get(h,'parent'),'color','w');



function HtSavePlot(figHandle, Options, suffix)
% function HtSavePlot(figHandle, Options, suffix)

if ~isempty(Options.PlotStyle.nameString)
    saveName = ['Track_',Options.PlotStyle.nameString,'_',suffix];
    printAll(figHandle,saveName);
end



function Ridge = HtRidgeFollow(EdgeSegments, Options, imgKev)
% function Ridge = HtRidgeFollow(EdgeSegments, Options, imgKev)
% 
% Follow the "ridge" of an electron track image.
% 
% Inputs: {standard choice}
%   EdgeSegments: structure requiring fields:
%       startRidgePix: coordinates of img to start at, in pixels
%       startDirectionIndices: angle to start moving in, units of 
%           Options.cutAngleIncrementDegrees
%   Options: structure requiring fields:
%       cutInterpolationMethod: for interp2. {'linear'}, 'spline'
%       trackEndLowThresholdKev: for ending the ridge-following. {0.1}
%       infiniteLoopThresholdPix: minimum distance from a previous point, 
%           in order to stop ridge-following because of infinite loop.
%           {0.5 * Options.xxxxxxx}
%       cutAngleIncrementDegrees: angular resolution of step direction.
%           (the actual recorded direction is after centroid adjustment.)
%       twoPiIndices: 360 degrees, in units of Options.cutAngleIncrementDegrees
%       searchAngleIndices: maximum turn per step = 0.5*searchAngleIndices.
%       cutCoordinatesPix: cell array of vectors. Each represents the
%           coordinates of one angled cut, relative to the ridge point.
%       cutDistanceFromCenter: the distance from center for cutCoordinatesPix
%          vectors.
%       pixelSizeUm: the pixel pitch of the image.
%       cutSamplingIntervalPix: the increment between neighboring
%           cutCoordinatesPix points.
%       cutAngleDegrees: [degrees] = cutAngleDegrees([indices])
%       positionStepSizePix: nominal distance from one ridge point to the next.
% 
% Outputs:
%   Ridge: structure containing the following fields:
%       positionPix: n-by-2 array containing coordinate pairs of ridge points.
%       fwhmUm: measured FWHM of the chosen cut, in um.
%       dedxKevUm: energy integrated over cut, into units of keV/um.
%       stepLengthPix: actual distance from the previous ridge point
%       alphaDegrees: actual angle from the previous ridge point

% set values for first step
Ridge.positionPix(1,1:2) = EdgeSegments.startRidgePix ;
startDirectionIndices = EdgeSegments.startDirectionIndices;
isFinished = false;

while ~isFinished
    Ridge = HtTakeOneStep(imgKev, Ridge, Options, startDirectionIndices);
    
    if Options.PlotStyle.single
        Options.PlotStyle.h = HtPlotRidgeStep(imgKev, Ridge, Options);
    end
    
    % are we at the end of the track?
    thisXPix = Ridge.positionPix(end,1);
    thisYPix = Ridge.positionPix(end,2);
    % there is a possibility (especially for noisy images)
    %   that all the cut points were under threshold, and
    %   this*Pix are NaN's. This happens at the end of the track.
    if ~isnan(thisXPix) && ~isnan(thisYPix)
        thisEnergyKev = interp2(imgKev, thisYPix, thisXPix, ...
                            Options.cutInterpolationMethod);
        isAtEndOfTrack = thisEnergyKev < Options.trackEndLowThresholdKev;
    else
        isAtEndOfTrack = true;
        lastGoodRidgePoint = find(~isnan(Ridge.positionPix(:,1)),1,'last');
        % erase NaN points
        Ridge.positionPix(lastGoodRidgePoint+2:end,:) = []; % 2 b/c trimmed later
        Ridge.fwhmUm(lastGoodRidgePoint+1:end) = [];
        Ridge.dedxKevUm(lastGoodRidgePoint+1:end) = [];
        Ridge.stepLengthPix(lastGoodRidgePoint+1:end) = [];
        Ridge.alphaDegrees(lastGoodRidgePoint+1:end) = [];
        Ridge.bestCutCoordinatesPix = ...
            Ridge.bestCutCoordinatesPix(1:lastGoodRidgePoint);  % cell array
    end
    % are we stuck in an infinite loop?
    % don't use the current point... it's going to move. use previous point.
    previousPoints = 1 : size(Ridge.positionPix,1)-2;
    prevXPix = Ridge.positionPix(end-1,1);
    prevYPix = Ridge.positionPix(end-1,2);
    dx = prevXPix - Ridge.positionPix(previousPoints,1);
    dy = prevYPix - Ridge.positionPix(previousPoints,2);
    threshold = Options.infiniteLoopThresholdPix;
	isInInfiniteLoop = any(dx.^2 + dy.^2 < threshold^2);
    
    % set values for next step
    isFinished = isAtEndOfTrack || isInInfiniteLoop;
    startDirectionIndices = []; %only use it for the first step
end

if isInInfiniteLoop
    Ridge.err = 'Infinite loop';
end
% trim last point
Ridge.positionPix(end,:) = [];

% plots
if Options.PlotStyle.multiple(11)
    % all ridge points
    [~, Options.PlotStyle.axesPosition] = ...
        HtPlotImage(imgKev, Options.PlotStyle.cmap, ...
        Options.PlotStyle.axesPosition, ...
        Options.PlotStyle.insetPix);
    hold on
    plot3(Ridge.positionPix(:,2)+0.5, Ridge.positionPix(:,1)+0.5, ...
        Options.PlotStyle.z0*ones(size(Ridge.positionPix(:,1))), ...
        'd','markerfacecolor',Options.PlotStyle.ptColor, ...
        'markeredgecolor',Options.PlotStyle.ptColor, ...
        'markersize',4);
    HtSavePlot(gcf, Options, '11_finished_points');
end
if Options.PlotStyle.multiple(12)
    % all ridge cuts
    [~, Options.PlotStyle.axesPosition] = ...
        HtPlotImage(imgKev, Options.PlotStyle.cmap, ...
        Options.PlotStyle.axesPosition, ...
        Options.PlotStyle.insetPix);
    hold on
    for i=1:length(Ridge.bestCutCoordinatesPix)
        z0 = Options.PlotStyle.z0 * ...
            ones(size(Ridge.bestCutCoordinatesPix{i}(:,1)));
        plot3(Ridge.bestCutCoordinatesPix{i}(:,2)+0.5, ...
            Ridge.bestCutCoordinatesPix{i}(:,1)+0.5, ...
            z0, '-', 'color',Options.PlotStyle.bestCutColor);
    end
    HtSavePlot(gcf, Options, '12_finished_cuts');
end




function Ridge = HtTakeOneStep(preparedImageKev, Ridge, Options, ...
    startDirectionIndices)
% function Ridge = HtTakeOneStep(preparedImageKev, Ridge, Options, ...
%     startDirectionIndices)
% 
% Take one step of ridge-following, making the accompanying measurements.
% 
% If it is the first step, startDirectionIndices must be supplied.
% 
% Ridge.position(end,:) should contain the point to start from and measure.
% That point will be modified to fit the centroid of the energy cut.
% Other Ridge properties should be one index shorter, because they have not been
%   measured yet.

thisRidgePointIndex = size(Ridge.positionPix,1);
previousRidgePointIndex = thisRidgePointIndex - 1;  %may or may not exist

thisXPix = Ridge.positionPix(thisRidgePointIndex,1);
thisYPix = Ridge.positionPix(thisRidgePointIndex,2);

if (nargin==3 || isempty(startDirectionIndices)) && previousRidgePointIndex > 0
    % startDirectionIndices not supplied
    % so, read it from the previous step
    lastStepPix1 = Ridge.positionPix(thisRidgePointIndex, 1:2) - ...
                   Ridge.positionPix(previousRidgePointIndex, 1:2);
    % lastStep will change once we adjust to centroid. hence lastStep1.
    startDirectionDegrees = 180/pi * atan2(lastStepPix1(2),lastStepPix1(1));
    startDirectionIndices = round(startDirectionDegrees / ...
        Options.cutAngleIncrementDegrees);
    % wrap around to positive angles
    if startDirectionIndices <= 0
        startDirectionIndices = startDirectionIndices + Options.twoPiIndices;
    end
elseif nargin ==2 && previousRidgePointIndex == 0
    % problem
    error('I need a startDirection!')
end
minimumCutAngleIndices = startDirectionIndices - Options.searchAngleIndices / 2;
maximumCutAngleIndices = startDirectionIndices + Options.searchAngleIndices / 2;
theseCutAnglesIndices = minimumCutAngleIndices:maximumCutAngleIndices;
%wrap around to positive angles in range
theseCutAnglesIndices(theseCutAnglesIndices <= 0) = ...
    theseCutAnglesIndices(theseCutAnglesIndices <=0) + ...
    Options.twoPiIndices;
theseCutAnglesIndices(theseCutAnglesIndices > Options.twoPiIndices) = ...
    theseCutAnglesIndices(theseCutAnglesIndices > Options.twoPiIndices) - ...
    Options.twoPiIndices;
%initialize
allCutsXPix                  = cell(1,length(theseCutAnglesIndices));
allCutsYPix                  = cell(1,length(theseCutAnglesIndices));
allCutsDistanceCoordinatePix = cell(1,length(theseCutAnglesIndices));
allCutsEnergyKev             = cell(1,length(theseCutAnglesIndices));
widthMetric                  =  nan(1,length(theseCutAnglesIndices));
for cutAngleNo = 1:length(theseCutAnglesIndices)
    thisCutXPix = thisXPix + ...
        Options.cutCoordinatesPix{theseCutAnglesIndices(cutAngleNo)}(:,1);
    thisCutYPix = thisYPix + ...
        Options.cutCoordinatesPix{theseCutAnglesIndices(cutAngleNo)}(:,2);
    
    % exclude out-of-bounds points
    minX = 1;
    minY = 1;
    maxX = size(preparedImageKev,1);
    maxY = size(preparedImageKev,2);
    isInBounds = thisCutXPix > minX & thisCutXPix < maxX & ...
                 thisCutYPix > minY & thisCutYPix < maxY;
    thisCutXPix = thisCutXPix(isInBounds);
    thisCutYPix = thisCutYPix(isInBounds);
    thisCutDistanceCoordinatePix = Options.cutDistanceCoordinatePix(isInBounds);
    thisCutDistanceFromCenterPix = Options.cutDistanceFromCenterPix(isInBounds);
    thisCutEnergyKev = interp2(preparedImageKev, thisCutYPix, thisCutXPix, ...
        Options.cutInterpolationMethod);
    
    % check for cutLowThreshold and exclude points
    pointIsIncluded = true(size(thisCutEnergyKev));
    % each side separately
    side1underThreshold = thisCutDistanceCoordinatePix < 0 & ...
            thisCutEnergyKev' < Options.cutLowThresholdKev;
    if any(side1underThreshold)
        endIndex = 1;
        pointsToExclude = endIndex:find(side1underThreshold,1,'last');
        pointIsIncluded(pointsToExclude) = false;
    end
    side2underThreshold = thisCutDistanceCoordinatePix > 0 & ...
            thisCutEnergyKev' < Options.cutLowThresholdKev;
    if any(side2underThreshold)
        endIndex = length(side2underThreshold);
        pointsToExclude = find(side2underThreshold,1,'first'):endIndex;
        pointIsIncluded(pointsToExclude) = false;
    end
    thisCutXPix = thisCutXPix(pointIsIncluded);
    thisCutYPix = thisCutYPix(pointIsIncluded);
    thisCutDistanceCoordinatePix = thisCutDistanceCoordinatePix(pointIsIncluded);
    thisCutDistanceFromCenterPix = thisCutDistanceFromCenterPix(pointIsIncluded);
    thisCutEnergyKev = thisCutEnergyKev(pointIsIncluded);
    
    % here, the previous algorithm excludes cut points less than 0.
    widthMetric(cutAngleNo) = sum(thisCutEnergyKev .* ...
        thisCutDistanceFromCenterPix');
    %save for centroid calculation
    allCutsXPix{cutAngleNo}                  = thisCutXPix;
    allCutsYPix{cutAngleNo}                  = thisCutYPix;
    allCutsDistanceCoordinatePix{cutAngleNo} = thisCutDistanceCoordinatePix;
    allCutsEnergyKev{cutAngleNo}             = thisCutEnergyKev;
end

[~,bestCutIndex] = min(widthMetric);
bestCutXPix                  = allCutsXPix{bestCutIndex};
bestCutYPix                  = allCutsYPix{bestCutIndex};
bestCutDistanceCoordinatePix = allCutsDistanceCoordinatePix{bestCutIndex};
bestCutEnergyKev             = allCutsEnergyKev{bestCutIndex};

% adjust to centroid
centroidXPix = sum(bestCutEnergyKev .* bestCutXPix) / sum(bestCutEnergyKev);
centroidYPix = sum(bestCutEnergyKev .* bestCutYPix) / sum(bestCutEnergyKev);
thisXPix = centroidXPix;
thisYPix = centroidYPix;

% measure stuff
thisFwhmUm = ...
    HtMeasureCutWidth(bestCutDistanceCoordinatePix * Options.pixelSizeUm, ...
    bestCutEnergyKev);
thisDedxKevUm = sum(bestCutEnergyKev) * Options.cutSamplingIntervalPix / ...
	Options.pixelSizeUm;
% thisAlpha measures from the previous point to here
if previousRidgePointIndex > 0
    % not the first point
    lastStepPix2 = Ridge.positionPix(thisRidgePointIndex, 1:2) - ...
                   Ridge.positionPix(previousRidgePointIndex, 1:2);
    thisAlphaDegrees = 180/pi * atan2(lastStepPix2(2), lastStepPix2(1));
    thisStepLengthPix = sqrt(lastStepPix2(1)^2 + lastStepPix2(2)^2);
else
    % first point only
    theseCutAnglesDegrees = Options.cutAngleDegrees(theseCutAnglesIndices);
    thisAlphaDegrees = theseCutAnglesDegrees(bestCutIndex);
    thisStepLengthPix = Options.positionStepSizePix;
    
    % plots
    if Options.PlotStyle.multiple(8)
        % image with all cuts through initial point
        [~, Options.PlotStyle.axesPosition] = ...
            HtPlotImage(preparedImageKev, Options.PlotStyle.cmap, ...
            Options.PlotStyle.axesPosition, ...
            Options.PlotStyle.insetPix);
        hold on;
        for i=1:length(theseCutAnglesDegrees)
            z0 = Options.PlotStyle.z0*ones(size(allCutsXPix{i}));
            plot3(allCutsYPix{i}+0.5, allCutsXPix{i}+0.5, z0, '-', ...
                'color',Options.PlotStyle.cutColor);
        end
        z0 = Options.PlotStyle.z0*ones(size(bestCutXPix));
        plot3(bestCutYPix+0.5, bestCutXPix+0.5, z0, '-', ...
            'color', Options.PlotStyle.bestCutColor);
        HtSavePlot(gcf, Options, '08_initialcuts_all');
    end
    if Options.PlotStyle.multiple(9)
        % image with three example cuts, sample points shown
        [~, Options.PlotStyle.axesPosition] = ...
            HtPlotImage(preparedImageKev, Options.PlotStyle.cmap, ...
            Options.PlotStyle.axesPosition, ...
            Options.PlotStyle.insetPix);
        hold on;
        colors{1}                   = Options.PlotStyle.cutColor;
        colors{bestCutIndex}        = Options.PlotStyle.bestCutColor;
        colors{length(allCutsXPix)} = Options.PlotStyle.otherCutColor;
        for i = [1, bestCutIndex, length(allCutsXPix)]
            z0 = Options.PlotStyle.z0*ones(size(allCutsXPix{i}));
            plot3(allCutsYPix{i}+0.5, allCutsXPix{i}+0.5, z0, 'd', ...
                'markerfacecolor', colors{i}, 'markersize', 4, ...
                'markeredgecolor', colors{i});
        end
        HtSavePlot(gcf, Options, '09_initialcuts_select');
    end
    if Options.PlotStyle.multiple(10)
        % graph of interpolated energy in each cut from (9)
        % first make axes object with proper settings...
        HtPlotImage(ones(2),'gray',Options.PlotStyle.axesPosition,[]);
        hold off
        colors{1}                   = Options.PlotStyle.cutColor;
        colors{bestCutIndex}        = Options.PlotStyle.bestCutColor;
        colors{length(allCutsXPix)} = Options.PlotStyle.otherCutColor;
        for i = [1, bestCutIndex, length(allCutsXPix)]
            plot(allCutsDistanceCoordinatePix{i}, allCutsEnergyKev{i}, ...
                '-o', 'color', colors{i}, ...
                'markeredgecolor','k', 'markerfacecolor', colors{i}, ...
                'linewidth',2, 'markersize',8);
            hold on;
        end
        xlabel('Length along cut [pixel lengths]')
        ylabel('Interpolated energy value [keV]')
        HtSavePlot(gcf, Options, '10_initialcuts_energy');
    end
end

% nextAlpha measures from here to where the next step will be
theseCutAnglesDegrees = Options.cutAngleDegrees(theseCutAnglesIndices);
nextAlphaDegrees = theseCutAnglesDegrees(bestCutIndex);
% take the step!
nextStepPix = Options.positionStepSizePix * ...
              [cosd(nextAlphaDegrees), sind(nextAlphaDegrees)];
nextXPix = thisXPix + nextStepPix(1);
nextYPix = thisYPix + nextStepPix(2);
nextRidgePointIndex = thisRidgePointIndex + 1;

% construct output structure
Ridge.positionPix          (thisRidgePointIndex,1:2) = [thisXPix, thisYPix];
Ridge.positionPix          (nextRidgePointIndex,1:2) = [nextXPix, nextYPix];
Ridge.fwhmUm               (thisRidgePointIndex)     = thisFwhmUm;
Ridge.dedxKevUm            (thisRidgePointIndex)     = thisDedxKevUm;
Ridge.stepLengthPix        (thisRidgePointIndex)     = thisStepLengthPix;
Ridge.alphaDegrees         (thisRidgePointIndex)     = thisAlphaDegrees;
Ridge.bestCutCoordinatesPix{thisRidgePointIndex}     = ...
    [bestCutXPix(:), bestCutYPix(:)];



function fwhm = HtMeasureCutWidth(cutXUm, cutEnergy)
% function fwhm = HtMeasureCutWidth(cutXUm, cutEnergy)
% 
% Return the full width at half maximum of an interpolated energy cut.
% 
% Inputs should be row vectors, with cutXUm in *microns*.

% Could either do a fit (traditional behavior), or a quick manual measurement.

% Check for empty vectors
if ~isempty(cutEnergy) && ~isempty(cutXUm)
    %{
    % Fit:
    minX = min(cutXUm);
    maxX = max(cutXUm);
    meanX = (minX + maxX) / 2;
    totalWidthX = maxX - minX;
    f = fit(cutXUm', cutEnergy', 'gauss1', ...
        'StartPoint', [max(cutEnergy),   meanX, totalWidthX/5], ...
        'Lower',      [max(cutEnergy)/2, minX,  0], ...
        'Upper',      [max(cutEnergy)*2, maxX,  totalWidthX]);
    fwhm = (f.c1 / sqrt(2)) * 2.355;
    %}
    % Manual: use fit_copy.m
    f = HtFitCopy(cutXUm', cutEnergy');
    fwhm = (f.c1 / sqrt(2)) * 2.355;
else
    fwhm = NaN;
end



function f = HtFitCopy(varargin)
%function f = HtFitCopy(xData, yData)
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



function h = HtPlotRidgeStep(imgKev, Ridge, Options)
% function h = HtPlotRidgeStep(imgKev, Ridge, Options)

OFFSET_PIX = 0.5;   % put the points in the center of the image pixels
Z = 1;              % vertical distance from image plane

if isfield(Options.PlotStyle,'h')
    h = Options.PlotStyle.h;
else
    %first step plotted for this track: make new axes with electron track image
    load cmaps
    h = SurfElectronTrack(imgKev,'cmap',cmaphotlog);
    Options.PlotStyle.h = h;
    hold on;
end

axes(Options.PlotStyle.h);  % bring to focus, make current

thisInd = length(Ridge.positionPix) - 1;
thisXY = Ridge.positionPix(thisInd,[2,1]) + ones(1,2)*OFFSET_PIX;
stepXY = Ridge.positionPix(thisInd+1,[2,1]) + ones(1,2)*OFFSET_PIX;
cutXY = Ridge.bestCutCoordinatesPix{thisInd}(:,[2,1]);
cutXY = cutXY + ones(size(cutXY))*OFFSET_PIX;

% ridge point
plot3(thisXY(1), thisXY(2), Z, '.g');
% cut
plot3(cutXY(:,1), cutXY(:,2), ones(size(cutXY,1),1)*Z, 'c');
% step direction
plot3([thisXY(1),stepXY(1)], [thisXY(2),stepXY(2)], [Z,Z], '-g');



function [Measurement,Ridge] = HtComputeDirection(trackEnergy, Ridge, Options)
% function [Measurement,Ridge] = HtComputeDirection(trackEnergy, Ridge, Options)

finalRidgePointIndex = size(Ridge.positionPix,1);   %already trimmed
reverseIndices = finalRidgePointIndex:-1:1;

Ridge.positionPix      = Ridge.     positionPix(reverseIndices,:);
Ridge.dedxKevUm        = Ridge.       dedxKevUm(reverseIndices);
Ridge.fwhmUm           = Ridge.          fwhmUm(reverseIndices);

% correction for alpha degrees
alpha_tmp = Ridge.alphaDegrees;
alpha_tmp(alpha_tmp < -135) = alpha_tmp(alpha_tmp < -135) + 360;
Ridge.alphaDegrees     = alpha_tmp(reverseIndices);
Ridge.stepLengthPix    = Ridge.   stepLengthPix(reverseIndices);

% 1. Measure width of track so we know how many points to skip.
preferredWidthMeasurementLengthPts = ...
    round(Options.widthMeasurementLengthPix / Options.positionStepSizePix);
if preferredWidthMeasurementLengthPts > 1 && ...
        preferredWidthMeasurementLengthPts <= finalRidgePointIndex
    actualWidthMeasurementLengthPts = preferredWidthMeasurementLengthPts;
elseif preferredWidthMeasurementLengthPts < 1
    actualWidthMeasurementLengthPts = 1;
elseif preferredWidthMeasurementLengthPts > finalRidgePointIndex
    actualWidthMeasurementLengthPts = finalRidgePointIndex;
end
measuredWidthUm = Options.measurementFunctionHandle(...
    Ridge.fwhmUm(1:actualWidthMeasurementLengthPts));
skipDiffusionPts = (measuredWidthUm - Options.pixelSizeUm) * ...
    4 / Options.pixelSizeUm;    % logbook 11/1/2009....
if skipDiffusionPts < 0
    skipDiffusionPts = 0;
end

cosBeta(1) = cosd(Options.initialBetaGuessDegrees);

% 2. get measurement selection range, assuming beta==0 and no diffusion
[parallelMeasurementStartPointNo, parallelMeasurementEndPointNo] = ...
    HtSelectMeasurementPoints(trackEnergy);

% 3. get measurement selection range, for a given beta, and skipping diffusion
actualMeasurementStartPointNo = ...
    ceil(parallelMeasurementStartPointNo * cosBeta(1) + skipDiffusionPts);
actualMeasurementEndPointNo = ...
    ceil(parallelMeasurementEndPointNo * cosBeta(1) + skipDiffusionPts);
if actualMeasurementStartPointNo > finalRidgePointIndex
    % fix start and end... this is kinda bad, only one measurement point
    actualMeasurementStartPointNo = finalRidgePointIndex;
    actualMeasurementEndPointNo = finalRidgePointIndex;
elseif actualMeasurementEndPointNo > finalRidgePointIndex
    actualMeasurementEndPointNo = finalRidgePointIndex;
end
if actualMeasurementEndPointNo < 1
    actualMeasurementStartPointNo = 1;
    actualMeasurementEndPointNo = 1;
elseif actualMeasurementStartPointNo < 1
    actualMeasurementStartPointNo = 1;
end
measurementIndices = actualMeasurementStartPointNo:actualMeasurementEndPointNo;

% 4. First estimate of beta, using selection which is calculated using beta0.
energyIndex = 1;
dedxIndex = 2;
if trackEnergy > Options.dedxTable(1,energyIndex) && ...
        trackEnergy < Options.dedxTable(end,energyIndex)
    dedxReference = interp1(Options.dedxTable(:,energyIndex), ...
                            Options.dedxTable(:,dedxIndex), ...
                            trackEnergy);
elseif trackEnergy <= Options.dedxTable(1,energyIndex)
    % energy too low to interpolate
    dedxReference = Options.dedxTable(1,dedxIndex);
elseif trackEnergy >= Options.dedxTable(end,energyIndex)
    % energy too high to interpolate
    dedxReference = Options.dedxTable(end,dedxIndex);
end
dedxMeasured = Options.measurementFunctionHandle(...
    Ridge.dedxKevUm(measurementIndices));
cosBeta(2) = dedxReference / dedxMeasured;
% closest physical solution
if cosBeta(2) > 1
    cosBeta(2) = 1;
end

% 5. Next selection calculation
actualMeasurementStartPointNo = ...
    ceil(parallelMeasurementStartPointNo * cosBeta(2) + skipDiffusionPts);
actualMeasurementEndPointNo = ...
    ceil(parallelMeasurementEndPointNo * cosBeta(2) + skipDiffusionPts);
if actualMeasurementStartPointNo > finalRidgePointIndex
    % fix start and end... this is kinda bad, only one measurement point
    actualMeasurementStartPointNo = finalRidgePointIndex;
    actualMeasurementEndPointNo = finalRidgePointIndex;
elseif actualMeasurementEndPointNo > finalRidgePointIndex
    actualMeasurementEndPointNo = finalRidgePointIndex;
end
if actualMeasurementEndPointNo < 1
    actualMeasurementStartPointNo = 1;
    actualMeasurementEndPointNo = 1;
elseif actualMeasurementStartPointNo < 1
    actualMeasurementStartPointNo = 1;
end
measurementIndices = actualMeasurementStartPointNo:actualMeasurementEndPointNo;

% 6. Second and final estimate of beta, using selection calculated from 
%       first estimate of beta.
dedxMeasured = Options.measurementFunctionHandle(...
    Ridge.dedxKevUm(measurementIndices));
cosBeta(3) = dedxReference / dedxMeasured;
% closest physical solution
if cosBeta(3) > 1
    cosBeta(3) = 1;
end

betaDegrees = acosd(cosBeta(3));

% Measure alpha using selection calculated from first estimate of beta
%       (same points as for final measurement of beta)

alphaDegrees = Options.measurementFunctionHandle(...
    Ridge.alphaDegrees(measurementIndices));
% That was the alpha pointing TOWARD the beginning of the track.
alphaDegrees = alphaDegrees + 180;
if alphaDegrees > 360
    alphaDegrees = alphaDegrees - 360;
end

% construct output
Measurement.alphaDegrees = alphaDegrees;
Measurement.betaDegrees = betaDegrees;
Measurement.dedxReference = dedxReference;
Measurement.dedxMeasured = dedxMeasured;
Measurement.indices = measurementIndices;
% Measurement.entrancePix = Ridge.positionPix(actualMeasurementStartPointNo, :);



function [measurementStartPointNo, measurementEndPointNo] = ...
    HtSelectMeasurementPoints(trackEnergy)
% function [measurementStartPointNo, measurementEndPointNo] = ...
%     HtSelectMeasurementPoints(trackEnergy)

% see node 1875 on bearing.berkeley.edu.
% "different indexing" refers to how the alpha value is saved in the ridge data.
%   old code had each step saving the angle toward the next step;
%   new code has each step saving the angle from the previous step.
% TODO: change this and see what happens.
measurementStartPointNo = sqrt(0.0825 * trackEnergy + 15.814) - 3.4;
measurementStartPointNo = measurementStartPointNo - 1;  % different indexing
measurementStartPointNo = max(measurementStartPointNo, 0);
measurementEndPointNo = measurementStartPointNo * 2 + 3.4;
measurementEndPointNo = measurementEndPointNo - 1;    % different indexing
measurementEndPointNo = max(measurementEndPointNo, 0);



function Output = HtSetOutput(preparedImageKev, Options, EdgeSegments, Ridge, ...
    Measurement)
% function Output = HtSetOutput(preparedImageKev, Options, EdgeSegments, Ridge, ...
%     Measurement)
Output.input = Options;
Output.img = preparedImageKev; % image with 0 padding
Output.energyT = Options.energyT;
Output.energyM = Options.energyM;
% Output.Etot = Options.energy; % E total, summing all the pixels up
% Output.ThSS = thetaStepSizeRadians;
Output.pixsize = Options.pixelSizeUm; % record the pixel size

Output.ends = length(EdgeSegments.energiesKev); % number of ends
Output.Eend = EdgeSegments.energiesKev(EdgeSegments.chosenIndex); % energy of the choosen ends
Output.thin = EdgeSegments.thinnedImage; % thinned image 
Output.lt = EdgeSegments.lowThresholdUsed; % threhold being used
Output.EdgeSegments = EdgeSegments; % ends structure

Output.x = Ridge.positionPix(:,1); % ridge row value
Output.y = Ridge.positionPix(:,2); % ridge col value
Output.w = Ridge.fwhmUm; % width of the chosen cut, for diffusion
Output.a0 = Ridge.alphaDegrees; % measured alpha along the ridge
Output.dE = Ridge.dedxKevUm; % measured dedx(s)
Output.Ridge = Ridge; % ridge structure

Output.alpha = Measurement.alphaDegrees; % estimated alpha 
Output.beta = Measurement.betaDegrees; % estimated beta
Output.Measurement = Measurement;   % measured structure
%dedxReference, dedxMeasured, indices




if ~isempty(Options.cheat)
    %pass the cheat structure back out.
    Output.cheat = Options.cheat;
    %already added to cheat structure:
    %cosbeta using actual electron energy (e.g. if detector captured full energy of electron)

    %this needs to be fixed. fields of 'cheat' have changed. and offsets must be included.
    % Track.cheat.enddist = sqrt((Track.x(1)-cheat.x)^2 + (Track.y(1)-cheat.y)^2);  %2D distance from first algorithm point to true electron start

end


