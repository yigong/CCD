function currentOut = Geant4TrackHandling7(varargin)
%function out = Geant4TrackHandling(datFilename,mode)
%function out = Geant4TrackHandling(matFilename,mode)
%function out = Geant4TrackHandling(g4matrix,mode)
%function out = Geant4TrackHandling(...,argID,arg)
%
% Modified 1/3/2013. Changed energy thresholds to simple 3-sigma.
%                    CCDsegment4 code is *not* implemented... too complicated.
%
% Geant4TrackHandling7.m:
%   changing structure so "cheat" output makes sense for multiple scattering.
%   lots of changes to "cheat".
%   Also, 'noiselist' argument turns into 'noise' which can be single or vector.
%   And 'pixelsize' argument can be single or a vector.
%
% Modified 9/12/2012. Adding New Mode of Handing Geant4 Tracking Data to DiffuseTrack6.m  : 
%       % New in DiffuseTrack6:
%       Restructuring inputs in the form of [pixel_plane_position, back_plane_position] rather 
%       than the assortment of flags used previously.
%
% Modified 5/1/2012. Added pixelsize argument to pass into DiffuseTrack5.
%
% Modified 4/13/2012. Fixed bug that skipped certain binding energy Edeps from multiple scattering of a single photon.
%                     Also cleaned up some oldcodeflag sections that could never be reached.
%                     Now known as Geant4TrackHandling5.m for compliance with Amy's numbering system.
% Modified 3/7/2012. May still have trouble loading dat or mat files with multiple events..
%
% Depends on:
%   DiffuseTrack6.m (or DiffuseTrack2b.m for diffusetrack2bflag)
%   AddCCDImageNoise.m
%   CCDsegment3.m
%   psft.mat, unless loaded as an optional argument
%
% datFilename: path/filename of a *.dat ASCII data file from geant4. (see format below)
% matFilename: path/filename of a *.mat Matlab data file with geant4 matrix in it. (see format below)
% g4matrix: numerical array of geant4 data. Should be 15 columns as follows:
%       Columns: (units)
%       1: Track ID
%       2: Parent ID
%       3: Step num
%       4: Charge                           (e-)
%       5,6,7: Initial position (x,y,z)     (mm)
%       8,9,10: Final position (x,y,z)      (mm)
%       11: Accumulated track length        (mm)
%       12: Step length                     (mm)
%       13: Final energy (after step)       (eV)
%       14: Initial energy (before step)    (eV)
%       15: Energy deposited (dE)           (eV)
%     default coordinate system is CCD depth along x dimension.
%
% mode: string to specify what type of simulation the g4matrix comes from.
%   'coinc': coincidence experiment
%           e.g.: coincBench_5mil_008.dat
%   'single': single CCD with photons incident upon it. (algorithm tests)
%           e.g.: eTracks_Si_1501.mat
% argID: string to identify optional argument
% arg: value of optional argument
%
% argID list:
%  psft: cell array of PSFs, or path/filename of psft.mat (defaults to psft.mat in current directory)
%  fullimage: set flag for saving the entire CCD image into the output structure 
%     (don't use in a loop, for memory's sake)
%  progressbar: set flag for using 'progressbar' function for multi-track operation
%  CCDsize: [x,y] size of the CCD in pixels. Default 726x1454 for coinc mode.
%     'var': set flag for adapting simulated size of device based on the event. 
%     default 'var' for 'single' mode.)
%  diffusetrack2b: set flag for using old diffusion code. (version 2b, which handles the secondary corrections)
%     *dropping support for this soon . . even now it may not work well . . 
%  noise: specify the sigma of the black level noise. *** in keV still, though fields are shown in eV
%     Can be a single number or a vector.
%     If empty or not specified, will use default measured noise level (see AddCCDImageNoise.m)
%  energywindow: supply [Emin,Emax] for events to process. units: keV (converted to eV for geant4 comarison)
%     The energy window applies to the initial energy of the primary electron 
%       of the source photon. Escaped energy is not subtracted.
%  driftdim: provide either {'x','y','z'} or {1,2,3} to define which dimension in geant4 is the thickness of the CCD.
%     default for coinc mode: z
%     default for single mode: x
%  pixelsize: specify the pixel pitch in um.
%     Can be a single number or a vector.
%     If empty or not specified, defaults to 10.5 um.
%  pixelnoisematrix: a logical array. 
%     if element (i,j) is true, then pixelSize(i) and noise(j) will be used together
%     if element (i,j) is false, then pixelSize(i) and noise(j) will not be used together
%  bypassPhotonErrorFlag: 
%     true to bypass "cannot find photon step for this primary electron" error message
% New in DiffuseTrack6:
%  driftdim:  driftDim set to {'x', 'y', 'z'} ==> sets indDepth in DiffuseTrack6.m to {1,2,3}
%  depthcoordinates: [pixel_plane_position, back_plane_position] -- Defines Detector 'Active Thickness'
%  % New handling of Geant4 tracking information for passing to Diffusion:
%  % layers: Location of [pixel_plane_position, back_plane_position], 
%                   with driftDim set to {'x', 'y', 'z'}   
%                   default [ -0.65 0 ]
% Still to come.....
% %%%%%%%%%%%%%%%%%% NEW Layer handeling ARG %%%%%%%%%%%%%%%%%%%%%%%%
% 
%
% OUTPUT:
%       out.T{struct} - Imgages -length(out.T)= # of Electron Track Images
%       out.E(vector) - CCD Image Energy - length(out.E) # of E-Track Images
%       out.cheat{struct} - Cheat Information for DiffuseTrack5.m 
%                           {1,1} Electron Trajectory info of 1st scattered electron
%                           (Needs to someday account for multi-CS events)
%                           (Currently Not using .cheat structure for
%                           analysis)
%           .cheat.Etot:    actual electron energy        (keV)
%           .cheat.Edep:    energy deposited in detector  (keV)
%           .cheat.x0:      coordinates of start of trajectory (um)
%           .cheat.x:       list of 3D coordinates that have been diffused. z=depth.
%           .cheat.dE:      list of energies diffused. corresponds to D.x.
%             ... and many more ...


%need to add:
%  thickness (and remove hardcoded references to 0.650)
%  something about PSFTs and diffusion


%% Input handling

if nargin<2
    error('Need input data and mode as input arguments')
end

%data
if ischar(varargin{1}) && ~isempty(strfind(varargin{1},'.dat'))
    g4Matrix = load(varargin{1},'-ascii');
elseif ischar(varargin{1}) && ~isempty(strfind(varargin{1},'.mat'))
    g4MatrixCurrent = load(varargin{1});
    if length(fieldnames(g4MatrixCurrent))==1
        tmp2 = fieldnames(g4MatrixCurrent);
        tmp = g4MatrixCurrent.(tmp2{1});
        g4Matrix = tmp{1};
        clear Mtemp tmp
    else
        error('don''t know how to handle multiple variables in data file')
    end
elseif isnumeric(varargin{1})
    g4Matrix = varargin{1};
end

%mode
if ~ischar(varargin{2})
    error('Mode should be a string')
elseif strcmpi(varargin{2},'coinc') || strcmpi(varargin{2},'coincidence')
    inputMode = 'coinc';
elseif strcmpi(varargin{2},'single') || strcmpi(varargin{2},'singleccd') || strcmpi(varargin{2},'single ccd')
    inputMode = 'single';
else
    error('Mode not recognized')
end

currentOut = [];
psft = [];
fullImageFlag = false;  %default
progressbarFlag = false;
oldCodeFlag = false;
noise = [];
energyWindow = [-Inf,Inf];
%{
%segmentation parameters . . this gets complicated for varying pixel sizes and noise.
segmentThresholdBase = 0.06;    %not an argument yet..
segmentThresholdSlope = (0.55 - segmentThresholdBase) / 10.5^2;   % keV/um^2 of pixelSize
%   set so that 10.5 um pixels give a standard 0.55 keV threshold
%}
pixelThresholdInverseSlope = 4 * 10.5^2;    
%   gives 4 for 10.5 um pixels using the following equation:
%     pixelThreshold = round(pixelThresholdInverseSlope / pixelSize^2);

segmentThresholdSigma = 3;  %3 sigma threshold, always.
pixelSize = 10.5; % um
thickness = 650;   %um
% default new input arguments:
driftDim = 3; % Cs137 (Spec & Coinc Sims)
depthCoordinates = [-0.65, 0];
bypassPhotonErrorFlag = false;


if strcmp(inputMode,'coinc')
    variableSizeFlag = false;
elseif strcmp(inputMode,'single')
    variableSizeFlag = true;
end
%optional arguments
i=3;
while true
    if i>nargin
        break
    elseif strcmpi(varargin{i},'psft')
        if nargin>i && ischar(varargin{i+1})
            tmp = load(varargin{i+1});
            psft = tmp.psft;
        elseif nargin>i && iscell(varargin{i+1})
            psft = varargin{i+1};
        elseif nargin==i
            warning('Psft input not found; looking in current directory')
        else
            error('Psft input must be cell array or filename')
        end
        i=i+2;
    elseif strcmpi(varargin{i},'fullimage') || strcmpi(varargin{i},'fullimg') || ...
            strcmpi(varargin{i},'fullimageflag') || strcmpi(varargin{i},'fullimgflag')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==1 || varargin{i+1}==0)))
            fullImageFlag = logical(varargin{i+1});
            i=i+2;
        else
            fullImageFlag = true;
            i=i+1;
        end
    elseif strcmpi(varargin{i},'nonoise') || strcmpi(varargin{i},'nonoiseflag')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==1 || varargin{i+1}==0)))
            noNoiseFlag = logical(varargin{i+1});
            i=i+2;
        else
            noNoiseFlag = true;
            i=i+1;
        end
        %not actually using this flag anymore . . use 'noise'
        if noNoiseFlag
            noise = 0;
        end
    elseif strcmpi(varargin{i},'progressbar') || strcmpi(varargin{i},'progressbarflag')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==1 || varargin{i+1}==0)))
            progressbarFlag = logical(varargin{i+1});
            i=i+2;
        else
            progressbarFlag = true;
            i=i+1;
        end
    elseif strcmpi(varargin{i},'variablesize') || strcmpi(varargin{i},'variableCCDsize')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==1 || varargin{i+1}==0)))
            variableSizeFlag = logical(varargin{i+1});
            i=i+2;
        else
            variableSizeFlag = true;
            i=i+1;
        end
    elseif strcmpi(varargin{i},'CCDsize')
        if nargin>i && isnumeric(varargin{i+1}) && length(varargin{i+1})==2
            CCDsize{1} = varargin{i+1};
            i=i+2;
        elseif nargin>i && ischar(varargin{i+1}) && (strcmpi(varargin{i+1},'var') || strcmpi(varargin{i+1},'variable'))
            variableSizeFlag = true;
            i=i+2;
        end
    elseif strcmpi(varargin{i},'pixelsize') || strcmpi(varargin{i},'pixsize')
        if nargin>i && isnumeric(varargin{i+1})
            pixelSize = varargin{i+1};
            i=i+2;
        end
    elseif strcmpi(varargin{i},'diffusetrack2b') || strcmpi(varargin{i},'oldcode')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==1 || varargin{i+1}==0)))
            oldCodeFlag = logical(varargin{i+1});
            i=i+2;
        else
            oldCodeFlag = true;
            i=i+1;
        end
    elseif strcmpi(varargin{i},'noise')
        if nargin>i && isnumeric(varargin{i+1})
            noise = varargin{i+1};
            i=i+2;
        elseif nargin>i && isempty(varargin{i+1})
            noise = [];
            i=i+2;
        else
            error('problem with noiselist input argument')
        end
    elseif strcmpi(varargin{i},'energywindow') || strcmpi(varargin{i},'energywindows')
        if nargin>i && isnumeric(varargin{i+1})
            energyWindow = varargin{i+1} * 1e3;
            i=i+2;
        else
            error('problem with energywindow input argument')
        end
    elseif strcmpi(varargin{i},'driftdim') || strcmpi(varargin{i},'depth')
        if nargin>i && isnumeric(varargin{i+1})
            driftDim = varargin{i+1};
            i=i+2;
        elseif nargin>i && ischar(varargin{i+1}) && length(varargin{i+1})==1
            if strcmpi(varargin{i+1},'x') || strcmpi(varargin{i+1},'1')
                driftDim = 1;
                i=i+2;
            elseif strcmpi(varargin{i+1},'y') || strcmpi(varargin{i+1},'2')
                driftDim = 2;
                i=i+2;
            elseif strcmpi(varargin{i+1},'z') || strcmpi(varargin{i+1},'3')
                driftDim = 3;
                i=i+2;
            else
                error('problem with driftdim input argument')
            end
        else
            error('problem with driftdim input argument')
        end
    elseif (strcmpi(varargin{i}, 'depthcoordinates') || ...
            strcmpi(varargin{i}, 'depthcoordinatesflag') || ...
            strcmpi(varargin{i}, 'depthcoord')); 
        if nargin>i && isnumeric(varargin{i+1}) && length(varargin{i+1})==2
            depthCoordinates = (varargin{i+1});
            i=i+2;
        else
            warning(['problem with depthcoordinates argument, ',...
                'should be 2 value vector [pixel_plane_position, back_plane_position]'])
            message = [ '... Using default depthcoordinates of: ' num2string(depthCoordinates)];
            warning(message);
            i=i+2;
        end
    elseif strcmpi(varargin{i}, 'pixelnoisematrix') || ...
            strcmpi(varargin{i}, 'pixelnoiselogical') || ...
            strcmpi(varargin{i}, 'pixelnoisearray')
        if nargin>i && (isnumeric(varargin{i+1}) || islogical(varargin{i+1}))
            pixelNoiseMatrix = varargin{i+1};
            i=i+2;
        else
            error('problem with pixelnoisematrix argument')
        end
        
        % Adding Layers ARG Option here: 9/28/2012 ABC
    elseif strcmpi(varargin{i}, 'bypassPhotonErrorFlag')
        if nargin>i && (isnumeric(varargin{i+1}) || islogical(varargin{i+1}))
            bypassPhotonErrorFlag = logical(varargin{i+1});
        else
            error('problem with bypassPhotonErrorFlag argument')
        end
        i=i+2;
    elseif (strcmpi(varargin{i}, 'layers') || ...
            strcmpi(varargin{i}, 'layersflag') || ...
            strcmpi(varargin{i}, 'layers')); 
        %
                
    else
        error(['Argument ID "',varargin{i},'" not recognized'])
        
    end
end

%% More initial setup.

%default psft
if isempty(psft)
    warning('PSF table not given; loading from file...')
    tmp = load('psft.mat');
    psft = tmp.psft;
end

%shortcuts to column numbers
    indTrackID = 1;
    indParentID = 2;
    indStepNum = 3;
    indCharge = 4;
    indInitPos = 5:7;
    indFinalPos = 8:10;
    indTrackLen = 11;
    indStepLen = 12;
    indFinalE = 13;
    indInitE = 14;
    inddE = 15;
    
%output structure
multiPixelSize = (length(pixelSize)>1); %logical
multiNoise = (length(noise)>1);         %logical

%default pixelnoisematrix
if isempty(who('pixelNoiseMatrix'))
    pixelNoiseMatrix = true(length(pixelSize),max(1,length(noise)));
end
%check pixelnoisematrix
if ~(size(pixelNoiseMatrix,1) == length(pixelSize) && ...
    size(pixelNoiseMatrix,2) == max(1,length(noise)))
    %problem
    error('pixelnoisematrix not the right size')
else
    pixelNoiseMatrix = logical(pixelNoiseMatrix);
end

%initialize
outputFieldsLength = sum(pixelNoiseMatrix(:));
outputFieldName = cell(1,outputFieldsLength);
fieldsPixelSizes = nan(1,outputFieldsLength);
fieldsNoises = nan(1,outputFieldsLength);

outputFieldNameIndex = 1;   %increment only according to pixelNoiseMatrix
if multiPixelSize && multiNoise
    %make cell array of field names, containing both pixel size value and noise value
    for pixelSizeIndex = 1:length(pixelSize)
        pixelSizeString = num2str(pixelSize(pixelSizeIndex));
        %replace the decimal point with an underscore, to make a valid field name
        if ~isempty(strfind(pixelSizeString,'.'))
            pixelSizeString(strfind(pixelSizeString,'.')) = '_';
        end
        
        for noiseIndex = 1:length(noise)
            if ~pixelNoiseMatrix(pixelSizeIndex,noiseIndex)
                %skip this combination
                continue
            end
            %show in eV
            noiseString = num2str(noise(noiseIndex).*1e3);
            if ~isempty(strfind(noiseString,'.'))
                noiseString(strfind(noiseString,'.')) = '_';
            end
            outputFieldName{outputFieldNameIndex} = ['pix',pixelSizeString,'noise',noiseString];
            %keep track of which fields have which pixel sizes
            fieldsPixelSizes(outputFieldNameIndex) = pixelSize(pixelSizeIndex);
            %keep track of which fields have which noise values
            fieldsNoises(outputFieldNameIndex) = noise(noiseIndex);
            outputFieldNameIndex = outputFieldNameIndex + 1;
        end
    end
elseif multiPixelSize && ~multiNoise
    %field names contain only pixel size value
    for pixelSizeIndex = 1:length(pixelSize)
        if ~pixelNoiseMatrix(pixelSizeIndex,1)
            %skip this one.. shouldn't encounter this case anyway
            continue
        end
        pixelSizeString = num2str(pixelSize(pixelSizeIndex));
        %replace the decimal point with an underscore, to make a valid field name
        if ~isempty(strfind(pixelSizeString,'.'))
            pixelSizeString(strfind(pixelSizeString,'.')) = '_';
        end
        outputFieldName{outputFieldNameIndex} = ['pix',pixelSizeString];
        %keep track of which fields have which pixel sizes
        fieldsPixelSizes(outputFieldNameIndex) = pixelSize(pixelSizeIndex);
        if isempty(noise)
            fieldsNoises(outputFieldNameIndex) = nan;   %this won't be seen later, anyway
        else
            fieldsNoises(outputFieldNameIndex) = noise; %one value of noise
        end
        outputFieldNameIndex = outputFieldNameIndex + 1;
    end
elseif ~multiPixelSize && multiNoise
    %field names contain only noise value
    for noiseIndex = 1:length(noise)
        if ~pixelNoiseMatrix(1,noiseIndex)
            %skip this one.. shouldn't encounter this case anyway
            continue
        end
        noiseString = num2str(noise(noiseIndex));
        if ~isempty(strfind(noiseString,'.'))
            noiseString(strfind(noiseString,'.')) = '_';
        end
        outputFieldName{outputFieldNameIndex} = ['noise',noiseString];
        %keep track of which fields have which noise values
        fieldsNoises(outputFieldNameIndex) = noise(noiseIndex);
        fieldsPixelSizes(outputFieldNameIndex) = pixelSize;
        outputFieldNameIndex = outputFieldNameIndex + 1;
    end
end
    
if isempty(who('CCDsize'))
%     CCDsize = [726,1454];   %pixels
    variableSizeFlag = true;
end
if variableSizeFlag
    CCDsize = cell(1,length(pixelSize));  %will be defined later...
end


% Mode: coincidence experiment simulation, OR single CCD.
%coordinates: center of CCD backplane is at (0,0,0); -z is into the detector.
if strcmp(inputMode,'coinc') || strcmp(inputMode,'single')
    %difference: coinc defaults to CCD depth in -z direction
    %            single defaults to CCD depth in +z direction (used to be +x direction)
    % *** remove 'mode' in a future version
    
    %first, check input arguments.
    
    %set indDepth. This is a number 1, 2 or 3 corresponding to x, y, or z.
    % indDepth is the dimension along the drift of the charge carriers, perpendicular to the pixel plane.
    if ~isempty(who('driftDim'))
        %specified in input arguments
        indDepth = driftDim;
    end
    
    %indLateral is whichver indices are not indDepth.
    indLateral = mod([indDepth+1,indDepth+2],3);
    indLateral(indLateral==0) = 3;
    % if indDepth = 3, then indLateral = [1,2];
    % if indDepth = 1, then indLateral = [2,3].

    %separate the source photons
    %   findNewSourcePhoton is the row index of the first step of each different source photon.
    findNewSourcePhoton = find(g4Matrix(:,indTrackID) == 1 ...
        & g4Matrix(:,indParentID) == 0 ...
        & g4Matrix(:,indStepNum) == 1 ...
        & g4Matrix(:,indCharge) == 0);
    %add an entry corresponding to the bottom row +1 of the matrix,
    %   to facilitate easy separation of the matrices
    findNewSourcePhoton(end+1) = size(g4Matrix,1)+1;
    
    %progressbar increments by source photon.
    if progressbarFlag
        progressbar(0);
    end
    
    %assign to output
    if length(findNewSourcePhoton)>2
        %multi-photon matrix: output will be a cell array of structures.
        %pre-allocate the cell array.
        out = cell(1,length(findNewSourcePhoton)-1);
%     else
%         %single photon matrix: output will be a structure. no preallocation needed
    end
    
    %% Each source photon.
    for sourcePhotonIndex = 1:length(findNewSourcePhoton)-1
        %get the section of matrix corresponding to this source photon (Event)
        g4MatrixCurrent = g4Matrix(findNewSourcePhoton(sourcePhotonIndex):...
            findNewSourcePhoton(sourcePhotonIndex+1)-1,:);
        
        if oldCodeFlag
            %use DiffuseTrack2b, which figures out how to do secondaries properly
            %   if secondary creation steps are not written
            D = DiffuseTrack2b(g4MatrixCurrent,psft,thickness,false);   
            %noiseflag in DiffuseTrack2b should be hardcoded to false
            
            %**oldCodeFlag section not set up for pixel size and noise value fieldnames
            currentCheat = D.cheat;
            if ~multiNoise && isempty(noise)
                %default to standard noise
                D.img = AddCCDImageNoise(D.img);
                [T,E] = CCDsegment3(D.img,segmentThreshold);
                if length(findNewSourcePhoton)==2
                    %single photon matrix: output a structure.
                    currentOut.T = T;
                    currentOut.E = E;
                else
                    %multi-photon matrix: output a cell array of structures.
                    currentOut.T = T;
                    currentOut.E = E;
                end
                currentOut.segmentThreshold = segmentThreshold;
            elseif ~multiNoise && ~isempty(noise)
                %single noise value, but not default
                D.img = AddCCDImageNoise(D.img,'blsig',noise(1));
                [T,E] = CCDsegment3(D.img,sqrt(segmentThreshold^2 + noise(1)^2));
                if length(findNewSourcePhoton)==2
                    %single photon matrix: output a structure.
                    currentOut.T = T;
                    currentOut.E = E;
                else
                    %multi-photon matrix: output a cell array of structures.
                    currentOut.T = T;
                    currentOut.E = E;
                end
                currentOut.segmentThreshold = segmentThreshold;
            else
                %multiple noise levels to simulate
                %must do everything multiple times
                currentOut = cell(1,length(noise));
                for noiseIndex=1:length(noise)
                    %add noise
                    imgtmp = AddCCDImageNoise(D.img,'blsig',noise(noiseIndex));
                    %segment
                    %use variable threshold: 0.55 and blsig, added in quadrature
                    threshold = sqrt(0.55^2 + noise(noiseIndex).^2);
                    [T,E] = CCDsegment3(imgtmp,threshold);
                    if length(findNewSourcePhoton)==2
                        %single photon matrix: one layer of cells, representing noise values
                        currentOut{noiseIndex}.T = T;
                        currentOut{noiseIndex}.E = E;
                        currentOut{noiseIndex}.cheat = currentCheat;
                        currentOut{noiseIndex}.blsig = noise(noiseIndex);
                        currentOut{noiseIndex}.th = threshold;
                    else
                        %multi-photon matrix: two layers of cells
                        currentOut{noiseIndex}.T = T;
                        currentOut{noiseIndex}.E = E;
                        currentOut{noiseIndex}.cheat = currentCheat;
                        currentOut{noiseIndex}.blsig = noise(noiseIndex);
                        currentOut{noiseIndex}.th = threshold;
                    end
                end
            end
            out{sourcePhotonIndex} = currentOut;
            continue
        end
        
        %define image bounds, for each pixel size.
        if variableSizeFlag
            %use extents of the geant4 track, with buffer around
            imageEdgeBufferMM = 210*1e-3;  %20 default pixels on each side. in mm.
            
            %ignore photons
            lgElectrons = g4MatrixCurrent(:,indCharge)==-1;
            
            %ignore escapes         %******is this redundant with below?
            lgContained = g4MatrixCurrent(:,indInitPos(indDepth)) > min(depthCoordinates) & ...
                g4MatrixCurrent(:,indInitPos(indDepth)) < max(depthCoordinates);
            
            lg = lgElectrons & lgContained;
            if ~any(lg)
                %no valid energy deposition here
                currentOut = [];
                continue
            end
            
            %in g4 units (mm).
            minX = min([g4MatrixCurrent(lg,indInitPos(indLateral(1))); g4MatrixCurrent(lg,indFinalPos(indLateral(1)))]);
            maxX = max([g4MatrixCurrent(lg,indInitPos(indLateral(1))); g4MatrixCurrent(lg,indFinalPos(indLateral(1)))]);
            minY = min([g4MatrixCurrent(lg,indInitPos(indLateral(2))); g4MatrixCurrent(lg,indFinalPos(indLateral(2)))]);
            maxY = max([g4MatrixCurrent(lg,indInitPos(indLateral(2))); g4MatrixCurrent(lg,indFinalPos(indLateral(2)))]);
            
            imageBoundsMM = [minX - imageEdgeBufferMM, maxX + imageEdgeBufferMM; ...
                minY - imageEdgeBufferMM, maxY + imageEdgeBufferMM];            %in mm
            
            %initialize
            imageBoundsPix = cell(1,length(pixelSize));
            CCDsize = cell(1,length(pixelSize));
            %each pixel size
            for pixelSizeIndex=1:length(pixelSize)
                imageBoundsPix{pixelSizeIndex} = round(imageBoundsMM *1e3./ pixelSize(pixelSizeIndex));  %in pixels
                
                %redefine CCD size.
                CCDsize{pixelSizeIndex} = (imageBoundsPix{pixelSizeIndex}(:,2) ...
                    - imageBoundsPix{pixelSizeIndex}(:,1) + [1;1])';
            end
        else
            %use real extents of CCD.
            %assume origin is at center of CCD in lateral dimensions.
            imageBoundsMM = 1e-3*[-CCDsize{1}(1)/2*pixelSize, CCDsize{1}(1)/2*pixelSize; ...
                -CCDsize{1}(2)/2*pixelSize, CCDsize{1}(2)/2*pixelSize];       %in mm
            %each pixel size
            for pixelSizeIndex=1:length(pixelSize)
                imageBoundsPix{pixelSizeIndex} = round(imageBoundsMM *1e3./ pixelSize(pixelSizeIndex));     %in pixels
            end
        end
        
        %define empty image for each pixel size.
        %initialize
        imageFinal = cell(1,length(pixelSize));
        for pixelSizeIndex=1:length(pixelSize)
            imageFinal{pixelSizeIndex} = zeros(CCDsize{pixelSizeIndex});
        end
        
        %identify primary electrons (daughters of source photon)
        sourcePhotonID = g4MatrixCurrent(find(g4MatrixCurrent(:,indCharge)==0,indTrackID),1);
        lgPrimaryElectron = g4MatrixCurrent(:,indParentID)==sourcePhotonID;
        primaryElectronList = unique(g4MatrixCurrent(lgPrimaryElectron,indTrackID));
        %unique also sorts particles into ascending order. cool!
        
        %initialize
        findXRay = nan(1,length(primaryElectronList));
        
        %% Each primary electron.
        for primaryElectronIndex = 1:length(primaryElectronList)
%             tic;
            %first step of this electron
            findFirstStep = find(g4MatrixCurrent(:,indTrackID)==primaryElectronList(primaryElectronIndex),1);
            firstStepVector = ...
                g4MatrixCurrent(findFirstStep,indFinalPos(1:3)) - ...
                g4MatrixCurrent(findFirstStep,indInitPos(1:3));
            %workaround for geant4 number precision issue: look at a longer step
            if size(g4MatrixCurrent,1)>findFirstStep ...
                    && g4MatrixCurrent(findFirstStep+1, indTrackID) == primaryElectronList(primaryElectronIndex) ...
                    && g4MatrixCurrent(findFirstStep+1, indStepNum) > 3
                %vector to beginning of step 10 (usually)
                %step number should be more than 3, to make sure it is far enough for precision not to matter.
                %   yeah this is kind of arbitrary.
                longStepVector = ...
                    g4MatrixCurrent(findFirstStep+1, indInitPos(1:3)) - ...
                    g4MatrixCurrent(findFirstStep, indInitPos(1:3));    %in mm
                longStepLength = g4MatrixCurrent(findFirstStep+1, indTrackLen) - ...
                    g4MatrixCurrent(findFirstStep+1, indStepLen) * 1e-3;   %um to mm
            elseif size(g4MatrixCurrent,1)>findFirstStep && ...
                    g4MatrixCurrent(findFirstStep+1, indTrackID) == primaryElectronList(primaryElectronIndex)
                %vector to end of that step
                %   I don't want to deal with finding subsequent steps that are further . . 
                %   this is a very small fraction of events to consider.
                longStepVector = ...
                    g4MatrixCurrent(findFirstStep+1, indInitPos(1:3)) - ...
                    g4MatrixCurrent(findFirstStep, indInitPos(1:3));    %mm
                longStepLength = g4MatrixCurrent(findFirstStep+1, indTrackLen) * 1e-3;  %um to mm
            else
                %single-step primary electron: not enough information to measure a longer step vector
                %   don't care about these anyway, this must be a very low-energy event
                longStepVector = nan(1,3);
                longStepLength = nan;
            end
            
            %initialize cheat structure.
            currentCheat.Edep = 0;     %so we can add each piece of energy deposition on, as it comes.
            currentCheat.dE = [];  %this gets concatenated with real info
            currentCheat.x0 = [];   %x0 is relative to true origin
            currentCheat.x = [];    %x is relative to true origin
            
            %known info
            currentCheat.particleID = primaryElectronList(primaryElectronIndex);
            currentCheat.E0 = g4MatrixCurrent(findFirstStep,indInitE) * 1e-3;   %eV to keV
            currentCheat.firstStepVector = firstStepVector;
            currentCheat.longStepLength = longStepLength;
            %define initial direction in CCD coordinates
            %   this assumes that DiffuseTrack does not rearrange lateral coordinates.
            %   In that case, alpha and beta should match the output of HybridTrack.
            %   But the surf/pcolor plots (including SurfElectronTrack) will be mirrored.
            %   And experimental tracks are flipped between x and y.
            currentCheat.alpha = atan2(firstStepVector(indLateral(2)), ...
                firstStepVector(indLateral(1))) * 180/pi;
            currentCheat.alphaLong = atan2(longStepVector(indLateral(2)), ...
                longStepVector(indLateral(1))) * 180/pi;
            %beta positive is away from the pixel plane
            if depthCoordinates(1) < depthCoordinates(2)
                currentCheat.beta = atan(firstStepVector(indDepth) / ...
                    sqrt(sum(firstStepVector(indLateral).^2))) * 180/pi;
                currentCheat.betaLong = atan(longStepVector(indDepth) / ...
                    sqrt(sum(longStepVector(indLateral).^2))) * 180/pi;
            else
                currentCheat.beta = atan(-firstStepVector(indDepth) / ...
                    sqrt(sum(firstStepVector(indLateral).^2))) * 180/pi;
                currentCheat.betaLong = atan(-longStepVector(indDepth) / ...
                    sqrt(sum(longStepVector(indLateral).^2))) * 180/pi;
            end
            
            %check energy window before going any further
            if currentCheat.E0 < energyWindow(1) || currentCheat.E0 > energyWindow(2)
                continue
            end
            
            %% Parent X-ray diffusion.
            %To identify the correct x-ray, we need to know the parent step.
            %   if multiplicity>1 (i.e. multiple Compton scatters), find the photon interaction
            %   closest to the start of this electron
            primaryElectronStartPosition = g4MatrixCurrent(findFirstStep,indInitPos);
            findSourcePhotonInteraction = find(g4MatrixCurrent(:,indTrackID)==sourcePhotonID);
            %   one would think that dE is always nonzero . . but not true. unsure why.
            
            for previousXRays = 1:length(findXRay)
                %if an x-ray was already used, it's not available
                findSourcePhotonInteraction(...
                    findSourcePhotonInteraction==findXRay(previousXRays)) = [];
            end
            photonInteractionDistanceSquared = ...
                sum((g4MatrixCurrent(findSourcePhotonInteraction,indFinalPos) - ...
                repmat(primaryElectronStartPosition,length(findSourcePhotonInteraction),1)).^2, 2);
            [minDistance, minDistanceIndex] = min(photonInteractionDistanceSquared);
            
%             disp(['xray distance is ',num2str(minDistance),' mm'])
            %for debug and verification
            currentCheat.XrayDistance = minDistance;
            
            if minDistance > 1e-3   %1 um...
                if bypassPhotonErrorFlag
                    %this just means that the photon interacted outside of the CCD, and then the electron scattered in.
                    currentCheat.Exray = nan;
                    currentCheat.x0 = primaryElectronStartPosition;
                    currentCheat.Etot = currentCheat.E0;   %already keV
                    currentCheat.sourcePhotonE1 = nan;
                    currentCheat.sourcePhotonE2 = nan;
                    currentCheat.sourcePhotonDirection1 = nan;
                    currentCheat.sourcePhotonDirection2 = nan;
                else
                    error('Cannot find the photon interaction for this primary electron');
                end
            else
                %save x-ray position for each primary electron, so we don't duplicate it.
                findXRay(primaryElectronIndex) = findSourcePhotonInteraction(minDistanceIndex);
                
                %Define photon information
                currentCheat.Exray = g4MatrixCurrent(findXRay(primaryElectronIndex),inddE) * 1e-3;  %eV to keV
                currentCheat.x0 = g4MatrixCurrent(findXRay(primaryElectronIndex),indFinalPos);
                currentCheat.Etot = currentCheat.E0 + currentCheat.Exray;   %already keV
                currentCheat.sourcePhotonE1 = g4MatrixCurrent(findXRay(primaryElectronIndex),indInitE) * 1e-3;  %eV to keV
                currentCheat.sourcePhotonE2 = g4MatrixCurrent(findXRay(primaryElectronIndex),indFinalE) * 1e-3; %eV to keV
                currentCheat.sourcePhotonDirection1 = ...
                    (g4MatrixCurrent(findXRay(primaryElectronIndex),indFinalPos) - ...
                    g4MatrixCurrent(findXRay(primaryElectronIndex),indInitPos)) ./ ...
                    sqrt(sum((g4MatrixCurrent(findXRay(primaryElectronIndex),indFinalPos) - ...
                    g4MatrixCurrent(findXRay(primaryElectronIndex),indInitPos)).^2));
                if g4MatrixCurrent(findXRay(primaryElectronIndex)+1,indTrackID)==sourcePhotonID
                    currentCheat.sourcePhotonDirection2 = ...
                        (g4MatrixCurrent(findXRay(primaryElectronIndex)+1,indFinalPos) - ...
                        g4MatrixCurrent(findXRay(primaryElectronIndex)+1,indInitPos)) ./ ...
                        sqrt(sum((g4MatrixCurrent(findXRay(primaryElectronIndex)+1,indFinalPos) - ...
                        g4MatrixCurrent(findXRay(primaryElectronIndex)+1,indInitPos)).^2));
    %                 cheat.sourcePhotonScatterAngle = 
                else
                    currentCheat.sourcePhotonDirection2 = nan(1,3);
                end
            
                %check active volume, then diffuse the x-ray.
                if g4MatrixCurrent(findXRay(primaryElectronIndex),indFinalPos(indDepth)) > min(depthCoordinates) && ...      zmin
                        g4MatrixCurrent(findXRay(primaryElectronIndex),indFinalPos(indDepth)) < max(depthCoordinates) && ...      zmax
                        g4MatrixCurrent(findXRay(primaryElectronIndex),indFinalPos(indLateral(1))) > imageBoundsMM(1,1) && ...    xmin
                        g4MatrixCurrent(findXRay(primaryElectronIndex),indFinalPos(indLateral(1))) < imageBoundsMM(1,2) && ...    xmax
                        g4MatrixCurrent(findXRay(primaryElectronIndex),indFinalPos(indLateral(2))) > imageBoundsMM(2,1) && ...    ymin
                        g4MatrixCurrent(findXRay(primaryElectronIndex),indFinalPos(indLateral(2))) < imageBoundsMM(2,2) && ...    ymax
                        g4MatrixCurrent(findXRay(primaryElectronIndex),inddE) > 0;     %nonzero deposition
                
                    %Diffuse the x-ray.
                    [imageFinal, currentCheat] = ...
                        DiffuseAndAdd(psft, g4MatrixCurrent(findXRay(primaryElectronIndex),:), ...
                        pixelSize, driftDim, depthCoordinates, multiPixelSize, imageFinal, imageBoundsPix, ...
                        currentCheat);
                end
    %             disp(['Diffusing photon xray took ',num2str(toc),' seconds'])
            end
            
%             tic;
            %% Each electron following from, and including, this primary electron.
            
            %get all particle IDs, then separate into electrons and brems.
            particleList = primaryElectronList(primaryElectronIndex);
            while true
                %hold on to the number of particles found...
                listLength = length(particleList);
                %find daughters of any particle in the list already
                for k = 1:listLength
                    %to qualify, particle must be a daughter
                    lgDaughters = g4MatrixCurrent(:,indParentID)==particleList(k);
                    %then add all of them to the list
                    particleList = [particleList, ...
                        unique(g4MatrixCurrent(lgDaughters,indTrackID))'];       %#ok<AGROW>
                end
                %remove all duplicate entries
                particleList = unique(particleList);
                %check if we're done: no new daughters left
                if listLength==length(particleList)
                    break
                end
                %otherwise: get another layer of daughters
            end
            
            %separate particleList into electrons and photons (brems)
            electronList = particleList;    %then we cut out brems
            bremsList = particleList;       %then we cut out electrons
            for particleIndex = 1:length(particleList)
                findFirstStep = find(g4MatrixCurrent(:,indTrackID)==particleList(particleIndex),1);
                if g4MatrixCurrent(findFirstStep,indCharge)==0
                    %this is a photon; mark for removal from electron list
                    electronList(particleIndex) = nan;
                else
                    %this is an electron; mark for removal from brems list
                    bremsList(particleIndex) = nan;
                end
            end
            %remove marked particles from lists
            electronList = electronList(~isnan(electronList));
            bremsList = bremsList(~isnan(bremsList));
            
            %for each electron
            for electronIndex = 1:length(electronList)
                %mark this particle's deposition that starts and ends within the active volume
                %   equal signs keep steps that end at the boundary.
                lgCCDactiveVolume = g4MatrixCurrent(:,indTrackID) == electronList(electronIndex) & ...  current particle
                    g4MatrixCurrent(:,indFinalPos(indDepth)) >= min(depthCoordinates) & ...                     zmin
                    g4MatrixCurrent(:,indFinalPos(indDepth)) <= max(depthCoordinates) & ...                     zmax
                    g4MatrixCurrent(:,indFinalPos(indLateral(1))) >= imageBoundsMM(1,1) & ...                   xmin
                    g4MatrixCurrent(:,indFinalPos(indLateral(1))) <= imageBoundsMM(1,2) & ...                   xmax
                    g4MatrixCurrent(:,indFinalPos(indLateral(2))) >= imageBoundsMM(2,1) & ...                   ymin
                    g4MatrixCurrent(:,indFinalPos(indLateral(2))) <= imageBoundsMM(2,2) & ...                   ymax
                    g4MatrixCurrent(:,indInitPos(indDepth)) >= min(depthCoordinates) & ...                      zmin
                    g4MatrixCurrent(:,indInitPos(indDepth)) <= max(depthCoordinates) & ...                      zmax
                    g4MatrixCurrent(:,indInitPos(indLateral(1))) >= imageBoundsMM(1,1) & ...                    xmin
                    g4MatrixCurrent(:,indInitPos(indLateral(1))) <= imageBoundsMM(1,2) & ...                    xmax
                    g4MatrixCurrent(:,indInitPos(indLateral(2))) >= imageBoundsMM(2,1) & ...                    ymin
                    g4MatrixCurrent(:,indInitPos(indLateral(2))) <= imageBoundsMM(2,2);    %                    ymax
                
                %check if we're done
                if ~any(lgCCDactiveVolume)
                    continue
                end
                
                %% Break this electron into segments as needed.
                
                g4MatrixElectron = g4MatrixCurrent(lgCCDactiveVolume,:);
                
                %hardcoded threshold for separating electron into separate segments. probably change for layers.
                dPositionsThreshold = 0.1;      %mm
                
                %check for large gaps in between electron segments: if electron escapes CCD and then re-enters
                %   this is to keep memory usage reasonable
                
                %first get 3-vector of position differential at each step [dx, dy, dz]
                dPositions = g4MatrixElectron(2:end,indInitPos) - g4MatrixElectron(1:end-1,indInitPos);
                %then get a distance metric: actually distance squared since we are just comparing it to a threshold.
                %   (square root is a relatively slow operation)
                dPositionsSquared = dPositions(:,1).^2 + dPositions(:,2).^2 + dPositions(:,3).^2;
                %where the distance metric is above the threshold, we are going to break the electron into segments.
                electronSegmentBreak = find(dPositionsSquared > dPositionsThreshold^2);
                
                %make separate matrices for each segment
                %   define cell array for each segment
                g4MatrixElectronSegment = cell(1,1+length(electronSegmentBreak));
                if ~isempty(electronSegmentBreak)
                    %split electron as determined above
                    g4MatrixElectronSegment{1} = g4MatrixElectron(1 : electronSegmentBreak(1), :);
                    for k=2:length(electronSegmentBreak)
                        g4MatrixElectronSegment{k} = g4MatrixElectron(electronSegmentBreak(k-1)+1 : electronSegmentBreak(k), :);
                    end
                    g4MatrixElectronSegment{end} = g4MatrixElectron(electronSegmentBreak(end)+1 : end, :);
                else
                    %no splitting needed
                    g4MatrixElectronSegment{1} = g4MatrixElectron;
                end
                
                %% Diffuse each segment.
                for segmentIndex = 1:length(g4MatrixElectronSegment)
                    [imageFinal, currentCheat] = ...
                        DiffuseAndAdd(psft, g4MatrixElectronSegment{segmentIndex}, ...
                        pixelSize, driftDim, depthCoordinates, multiPixelSize, imageFinal, imageBoundsPix, ...
                        currentCheat);
                end
            end
%             disp(['Diffusing all electrons took ',num2str(toc),' seconds'])
            
%             tic;
            %% Diffuse brems x-rays.
            currentCheat.Ebrems = 0;
            
            for bremsIndex = 1:length(bremsList)
                %add to Ebrems info, whether it is contained or not
                findFirstStep = find(g4MatrixCurrent(:,indTrackID)==bremsList(bremsIndex),1);
                currentCheat.Ebrems = currentCheat.Ebrems + g4MatrixCurrent(findFirstStep,indInitE) * 1e-3; %eV to keV
                
                %check active volume
                lgCCDactiveVolume = g4MatrixCurrent(:,indTrackID) == bremsList(bremsIndex) & ...       current particle
                    g4MatrixCurrent(:,indFinalPos(indDepth)) > min(depthCoordinates) & ...      zmin
                    g4MatrixCurrent(:,indFinalPos(indDepth)) < max(depthCoordinates) & ...      zmax
                    g4MatrixCurrent(:,indFinalPos(indLateral(1))) > imageBoundsMM(1,1) & ...    xmin
                    g4MatrixCurrent(:,indFinalPos(indLateral(1))) < imageBoundsMM(1,2) & ...    xmax
                    g4MatrixCurrent(:,indFinalPos(indLateral(2))) > imageBoundsMM(2,1) & ...    ymin
                    g4MatrixCurrent(:,indFinalPos(indLateral(2))) < imageBoundsMM(2,2) & ...    ymax
                    g4MatrixCurrent(:,inddE) > 0;                                   % deposited energy
                findCCDactiveVolume = find(lgCCDactiveVolume);
                
                %diffuse each interaction separately, if multiple scatters
                for bremsInteractionIndex = 1:length(findCCDactiveVolume)
                    [imageFinal, currentCheat] = ...
                        DiffuseAndAdd(psft, g4MatrixCurrent(findCCDactiveVolume(bremsInteractionIndex),:), ...
                        pixelSize, driftDim, depthCoordinates, multiPixelSize, imageFinal, imageBoundsPix, ...
                        currentCheat);
                end
            end
%             disp(['Diffusing all brems took ',num2str(toc),' seconds'])
            
            %save cheat structure for this primary electron
            %   do this for each configuration
            if ~multiNoise && ~multiPixelSize
                %only one configuration; no additional structure needed.
                currentOut.cheat(primaryElectronIndex) = currentCheat;
            else
                %each pixel size has a different cheat structure
                %for each pixel size
                for pixelSizeIndex = 1:length(pixelSize)
                    %see which configurations include this pixel size
                    findFieldNames = find(fieldsPixelSizes==pixelSize(pixelSizeIndex));
                    %for each configuration with this pixel size
                    for findFieldNamesIndex = 1:length(findFieldNames)
                        %assign the cheat structure associated with this pixel size
                        currentOut.(outputFieldName{findFieldNames(findFieldNamesIndex)}).cheat(primaryElectronIndex) = ...
                            currentCheat;
                    end
                end
            end
            
            %next primary electron
        end
        
        %final image noise and segmentation
        %configuration-dependent
        if ~multiNoise && ~multiPixelSize
            %simplest structure
            if isempty(noise)
                %default noise
                imageNoisy = AddCCDImageNoise(imageFinal{1});
                %add info
                currentOut.noise = 'default';
            else
                %specified noise
                imageNoisy = AddCCDImageNoise(imageFinal{1},'blsig',noise);
                %add info
                currentOut.noise = noise;
            end
            %add info
            currentOut.pixsize = pixelSize;
            %adjust threshold for pixel size. default threshold is for 10.5 um pixels.
            segmentThresholdCurrent = segmentThresholdSigma * noise;
            pixelThresholdCurrent = round(pixelThresholdInverseSlope / pixelSize^2);
            %segment or provide full image
            if fullImageFlag
                %no segmentation, just give the full image.
                %but with variable CCD size, this will still be a subimage.
                currentOut.img = imageNoisy;
            else
                %segment
                [currentTracks, currentOut.E] = CCDsegment3(imageNoisy, segmentThresholdCurrent, ...
                    'pixelthreshold', pixelThresholdCurrent);
                %adjust x,y
                for trackIndex = 1:length(currentTracks)
                    currentTracks{trackIndex}.x = currentTracks{trackIndex}.x + imageBoundsPix{1}(1);
                    currentTracks{trackIndex}.y = currentTracks{trackIndex}.y + imageBoundsPix{1}(2);
                end
                currentOut.T = currentTracks;
                
                currentOut.segmentThreshold = segmentThresholdCurrent;
                currentOut.pixelThreshold = pixelThresholdCurrent;
            end
        else
            %use field names defined up top
            
            for pixelSizeIndex = 1:length(pixelSize)
                if isempty(noise)
                    %default noise
                    imageNoisy = AddCCDImageNoise(imageFinal{pixelSizeIndex});
                    %identify which field we are writing to here
                    outputFieldNameIndex = find(fieldsPixelSizes==pixelSize(pixelSizeIndex));
                    if isempty(outputFieldNameIndex)
                        %disabled in pixelNoiseMatrix
                        continue
                    end
                    %set segmentation threshold for pixel size
                    error('I don''t know what the noise level is. AddCCDImageNoise won''t tell me.')
                    
                    pixelThresholdCurrent = round(pixelThresholdInverseSlope / pixelSize(pixelSizeIndex)^2);
                    %segment
                    if fullImageFlag
                        %no segmentation, just give the full image.
                        %but with variable CCD size, this will still be a subimage.
                        currentOut.(outputFieldName{outputFieldNameIndex}).img = imageNoisy;
                    else
                        %segment
                        [currentTracks, currentOut.(outputFieldName{outputFieldNameIndex}).E] = ...
                            CCDsegment3(imageNoisy, segmentThresholdCurrent, 'pixelthreshold', pixelThresholdCurrent);
                        %adjust x,y
                        for trackIndex = 1:length(currentTracks)
                            currentTracks{trackIndex}.x = currentTracks{trackIndex}.x + imageBoundsPix{pixelSizeIndex}(1);
                            currentTracks{trackIndex}.y = currentTracks{trackIndex}.y + imageBoundsPix{pixelSizeIndex}(2);
                        end
                        currentOut.(outputFieldName{outputFieldNameIndex}).T = currentTracks;
                        
                        currentOut.(outputFieldName{outputFieldNameIndex}).segmentThreshold = segmentThresholdCurrent;
                        currentOut.(outputFieldName{outputFieldNameIndex}).pixelThreshold = pixelThresholdCurrent;
                    end
                    %additional info
                    currentOut.(outputFieldName{outputFieldNameIndex}).pixsize = pixelSize(pixelSizeIndex);
                    %   pixelSize is saved in the DiffuseTrack7 output, 
                    %   but it isn't transferred into the TrackHandling data structure
                    currentOut.(outputFieldName{outputFieldNameIndex}).noise = 'default';
                    %   it would be nice to be able to write a value here, 
                    %   but that would require asking AddCCDImageNoise, or hardcoding it in here.
                    continue
                end
                %otherwise, we can loop through noise (even if there is only one value)
                for noiseIndex = 1:length(noise)
                    %get a noisy image from the right pixel size
                    imageNoisy = AddCCDImageNoise(imageFinal{pixelSizeIndex},'blsig',noise(noiseIndex));
                    %identify which field we are writing to here
                    outputFieldNameIndex = find(fieldsPixelSizes==pixelSize(pixelSizeIndex) ...
                        & fieldsNoises==noise(noiseIndex));
                    if isempty(outputFieldNameIndex)
                        %disabled in pixelNoiseMatrix
                        continue
                    end
                    %set segmentation threshold
                    if noise(noiseIndex)==0
                        %no noise: collect all the track (0.1 eV threshold)
                        segmentThresholdCurrent = 1e-4;
                    else
                        %noise: n sigma
                        segmentThresholdCurrent = segmentThresholdSigma * noise(noiseIndex);
                    end
                    %set pixel threshold
                    pixelThresholdCurrent = round(pixelThresholdInverseSlope / pixelSize(pixelSizeIndex)^2);
                    %segment
                    if fullImageFlag
                        %no segmentation, just give the full image.
                        %but with variable CCD size, this will still be a subimage.
                        currentOut.(outputFieldName{outputFieldNameIndex}).img = imageNoisy;
                    else
                        %segment
                        [currentTracks, currentOut.(outputFieldName{outputFieldNameIndex}).E] = ...
                            CCDsegment3(imageNoisy, segmentThresholdCurrent, 'pixelthreshold', pixelThresholdCurrent);
                        %adjust x,y
                        for trackIndex = 1:length(currentTracks)
                            currentTracks{trackIndex}.x = currentTracks{trackIndex}.x + imageBoundsPix{pixelSizeIndex}(1);
                            currentTracks{trackIndex}.y = currentTracks{trackIndex}.y + imageBoundsPix{pixelSizeIndex}(2);
                        end
                        currentOut.(outputFieldName{outputFieldNameIndex}).T = currentTracks;
                        
                        currentOut.(outputFieldName{outputFieldNameIndex}).segmentThreshold = segmentThresholdCurrent;
                        currentOut.(outputFieldName{outputFieldNameIndex}).pixelThreshold = pixelThresholdCurrent;
                    end
                    %additional info
                    currentOut.(outputFieldName{outputFieldNameIndex}).pixsize = pixelSize(pixelSizeIndex);
                    %   pixelSize is saved in the DiffuseTrack7 output, 
                    %   but it isn't transferred into the TrackHandling data structure
                    currentOut.(outputFieldName{outputFieldNameIndex}).noise = noise(noiseIndex);
                end
            end
        end
        
        if progressbarFlag
            %update progressbar
            progressbar(sourcePhotonIndex/(length(findNewSourcePhoton)-1));
        end
        
        %assign currentOut to out, depending on multiple source photons
        if length(findNewSourcePhoton)-1 > 1
            %cell array for source photons
            out{sourcePhotonIndex} = currentOut;
        else
            %no cell array
            out = currentOut;
        end
        
        %next source photon
    end
    
    if progressbarFlag
        %make sure progressbar is closed
        progressbar(1);
    end
    
end




function img = AddToImage(img,imageboundspix,D)
    %truncate at edge of CCD active area
    trunc(1) = max(0, imageboundspix(1,1) - D.offsets(1));  %left
    trunc(2) = max(0,-imageboundspix(1,2) + D.offsets(1) + size(D.img,1)); %right
    trunc(3) = max(0, imageboundspix(2,1) - D.offsets(2));  %bottom
    trunc(4) = max(0,-imageboundspix(2,2) + D.offsets(2) + size(D.img,2)); %top
    
    x1 = -imageboundspix(1,1) + D.offsets(1) + (trunc(1)+1:(size(D.img,1)-trunc(2)));
    y1 = -imageboundspix(2,1) + D.offsets(2) + (trunc(3)+1:(size(D.img,2)-trunc(4)));
    %add to main image
    img(x1, y1) = img(x1,y1) + D.img(trunc(1)+1:end-trunc(2), trunc(3)+1:end-trunc(4));
    

function [imageFinal, currentCheat] = DiffuseAndAdd(psft, g4MatrixSegment, ...
        pixelSize, driftDim, depthCoordinates, multiPixelSize, imageFinal, imageBoundsPix, ...
        currentCheat)
    
    %diffuse
    D = DiffuseTrack7(psft, 'geant4', g4MatrixSegment, ...
        'driftdim', driftDim, 'depthCoordinates', depthCoordinates, 'pixelSize', pixelSize);
    %output varies with number of pixel sizes
    if ~multiPixelSize
        %only one pixel size means a single structure from D

        %add into final image
        imageFinal{1} = AddToImage(imageFinal{1},imageBoundsPix{1},D);
        %extend cheat
        currentCheat.Edep = currentCheat.Edep + D.cheat.Edep;
        currentCheat.dE = [currentCheat.dE; D.cheat.dE];
        currentCheat.x = [currentCheat.x; D.cheat.x];
    else
        %multiple pixel sizes means a cell array from D
        for pixelSizeIndex = 1:length(pixelSize)
            %add into final image
            imageFinal{pixelSizeIndex} = AddToImage(imageFinal{pixelSizeIndex}, ...
                imageBoundsPix{pixelSizeIndex}, D{pixelSizeIndex});
        end
        %include x-ray in deposition information
        currentCheat.Edep = currentCheat.Edep + D{1}.cheat.Edep;
        currentCheat.dE = [currentCheat.dE; D{1}.cheat.dE];
        currentCheat.x = [currentCheat.x; D{1}.cheat.x];
    end