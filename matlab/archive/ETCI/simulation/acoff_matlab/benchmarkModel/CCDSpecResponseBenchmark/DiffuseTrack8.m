function [D] = DiffuseTrack8(varargin)
%function [D] = DiffuseTrack8(psft, mode, inputdata)
%function [D] = DiffuseTrack8(psft, mode, inputdata, argname, argvalue, ...)
%function [D] = DiffuseTrack8(psft, mode, inputdata, argflagname, ...)
%
% New in DiffuseTrack8:
%  Most of the 'cheat' structure is accomplished in TrackHandling,
%   so it is removed here.
%  Allows for multiple pixel sizes to be diffused in one function call to save computing time.
%   But they all must be multiples of psftLateralResolution.
%  Alpha and beta are calculated in Geant4TrackHandling now, not here.
%
% New in DiffuseTrack6:
%  Restructuring inputs in the form of [pixel_plane_position, back_plane_position] rather 
%   than the assortment of flags used previously.
%
% New in DiffuseTrack5:
%  Adding pixelsize as an input argument, and taking out all hardcoded numbers based on
%    10.5 um.
%
% New in DiffuseTrack4:
%  Designed to run within another script that sorts out particles and different layers.
%  The M input should ONLY include a single particle! not even secondaries.
%  Variable input arguments, including manual mode.
%
% psft: table of point spread functions, from MakePsfTable(tdt,D,0.5,zres,savename)
%        this should be a cell array of 121x121 psf's.
% mode: either "geant4" or "manual" (case-insensitive)
% inputdata (geant4 mode):
%       electron data matrix from geant4 for a SINGLE particle: M = g4M(g4M(:,1)==k,:)
%       Coordinate system: x_geant is depth dimension (pixel face = 0, back face = (th/1000))
%       (x,y) = (0,0) is the corner between 4 pixels, not the center of a pixel.
%       (coordinates need to be adjusted in parent code for detector systems with multiple CCDs)
%
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
% inputdata (manual mode):
%       list of (x,y,z,E) coordinates. 
%       (x,y,z) are in um.
%       z is depth dimension; z = 0 is pixel plane; max(z) should not exceed thickness.
%       (x,y) = (0,0) is the corner between 4 pixels, not the center of a pixel.
%       E is in keV.
%
% Optional arguments:
%   driftdim: dimension of charge drifting, either 'x' or 'y' or 'z' or 1 or 2 or 3.
%       defaults to 'z' (coincidence experiments, 9/2012 Cs-137 singles)
%   depthcoordinates: [pixel_plane_position, back_plane_position] in geant4 matrix units (mm)
%       defaults to [-0.65,0] (coincidence experiment).
%   pixelsize: pixel pitch, in um. (default 10.5)
%       this can be a list of pixel sizes to diffuse to.
%   energyaccounting: for use with Geant4TrackHandling, manual energy correction mode. 
%       uses InitE - FinalE, instead of dE, for intra-step deposition.
%
%output variables
%
%If pixel size is singular, or default:
%   D is a structure:
%   D.img:     data matrix of ccd plane pixels.    (keV)
%   D.offsets: absolute position of pixel (i,j) is (offsets(1)+i,offsets(2)+j). origin matches that of input data. 
%   D.pixsize: pixel size (um)
%   D.cheat:
%     cheat.Edep:    energy deposited in detector  (keV)
%     cheat.x:       list of 3D coordinates that have been diffused. z=depth. (um)
%     cheat.dE:      list of energies diffused. corresponds to D.x.
%
%If multiple pixel sizes are given:
%   D is a cell array of structures, each as above.
%
% future work: 
% change format of psft to include thickness or depth resolution.
%   currently it will be hardcoded to 650um thickness, so that we can apply it to 
%   a "thinner" but actually non-depleted CCD.
% make the fine grid resolution an optional argument (also depends on psft)

%% Input handling

%number of arguments
if nargin < 3
    error('Not enough input arguments: need psft, mode, inputdata')
end
%check first arg (psft)
if iscell(varargin{1})
    psft = varargin{1};
else
    error('psft should be a cell array of 121x121 matrices')
end
%check second arg (mode)
if ischar(varargin{2}) && strcmpi(varargin{2},'geant') || strcmpi(varargin{2},'geant4') ||...
        strcmpi(varargin{2},'g4') || strcmpi(varargin{2},'g')
    g4mode = true;
elseif ischar(varargin{2}) && strcmpi(varargin{2},'manual') || strcmpi(varargin{2},'m')
    g4mode = false;
else
    error('mode should be a string, either "geant4" or "manual"')
end
%check third arg (inputdata)
if isnumeric(varargin{3})
    inputMatrix = varargin{3};
else
    error('inputdata should be a numeric array')
end

%defaults
pixelSize = 10.5;
indDepth = 3;       %default to z dimension (coincidence geometry)
indLateral = mod([indDepth+1,indDepth+2],3);
indLateral(indLateral==0) = 3;
if g4mode
    depthCoordinates = [-0.65,0];   %default to coincidence geometry
else
    depthCoordinates = [0,0.65];    %default to positive
end
energyAccountingFlag = false;

%to be changed in future versions: get this info from the psft directly
psftThickness = 650;   %um
psftLateralResolution = 0.5; %um  %must match the resolution of the fine grid in PSFT
%um depth resolution of PSFT
psftDepthResolution = psftThickness/(length(psft)-1);
%I'm intentionally not scaling the depthcoordinates thickness to the 
%   psft_thickness . . . that would be bad.

%optional arguments
i=4;
while true
    if i > nargin
        break
    end
    if strcmpi(varargin{i},'pixelsize') || strcmpi(varargin{i},'pixelpitch') ...
            || strcmpi(varargin{i},'pixel') || strcmpi(varargin{i},'pixsize')
        if length(varargin)>i
            if isnumeric(varargin{i+1})
                pixelSize = varargin{i+1};
                i=i+2;
            else
                error('pixelsize should be numeric')
            end
        else
            warning('optional argument "pixelsize" not defined')
            i=i+1;
        end
    elseif strcmpi(varargin{i},'driftdim') || strcmpi(varargin{i},'drift') || strcmpi(varargin{i},'inddepth')
        if length(varargin)<=i
            error('optional argument "driftdim" not defined')
        end
        if isnumeric(varargin{i+1}) && length(varargin{i+1})==1 ...
                && varargin{i+1}>0 && varargin{i+1}<4
            indDepth = varargin{i+1};
        elseif ischar(varargin{i+1}) && length(varargin{i+1})==1 ...
                && (strcmpi(varargin{i+1},'x') || strcmpi(varargin{i+1},'y') ...
                || strcmpi(varargin{i+1},'z'))
            switch lower(varargin{i+1})
                case 'x'
                    indDepth = 1;
                case 'y'
                    indDepth = 2;
                case 'z'
                    indDepth = 3;
            end
        else
            error('driftdim input not recognized')
        end
        indLateral = mod([indDepth+1,indDepth+2],3);
        indLateral(indLateral==0) = 3;
        i=i+2;
    elseif strcmpi(varargin{i},'depthcoord') || ...
            strcmpi(varargin{i},'depthcoordinates') || ...
            strcmpi(varargin{i},'depthcoordinate') || ...
            strcmpi(varargin{i},'layers')
        if nargin==i
            error('Argument depthcoordinates not specified')
        end
        if ~isnumeric(varargin{i+1}) || length(varargin{i+1})-2~=0
            error('depthcoordinates should be a two-element vector')
        end
        depthCoordinates = varargin{i+1};   %mm
        i=i+2;
    elseif strcmpi(varargin{i},'energyaccountingflag') || ...
            strcmpi(varargin{i},'energyaccounting') || ...
            strcmpi(varargin{i},'manualenergycorrect')
        if nargin>i && (islogical(varargin{i+1}) || ...
                (isnumeric(varargin{i+1}) && (varargin{i+1}==0 || varargin{i+1}==1)))
            energyAccountingFlag = logical(varargin{i+1});
            i=i+2;
        else
            energyAccountingFlag = true;
            i=i+1;
        end
    else
        error(['Argument ID ',varargin{i},' not valid'])
    end
end

%% Fix up parameters.

%Convert depthcoordinates from mm to um
depthCoordinates = depthCoordinates*1e3;

%hold on to maximum pixel size
maxPixelSize = max(pixelSize);
if any(round(pixelSize/psftLateralResolution)~=(pixelSize/psftLateralResolution))
    %problem
    error('Pixel sizes not compatible with psft lateral resolution')
end

psfPatchWidth = size(psft{1},1);    %dimensions of psf, in units of the fine grid
psfPatchRadius = (psfPatchWidth-1)/2;

%Check for geometry errors
if (abs(depthCoordinates(2) - depthCoordinates(1)) - psftThickness) > 0    %allow rounding error?
    error('CCD geometry should not be thicker than CCD thickness described by diffusion table.')
end

%space around the energy deposition points, to accomodate plenty of diffusion.
%   this buffer is enough for any pixel size ~50 um or less, for standard PSFT.
imageEdgeBuffer = 85;    %um            %point (0,0) is in the center of a pixel

%% Set up energy deposition points.

if g4mode
    %shortcuts to column numbers
    indTrackID = 1;
%     indParentID = 2;
%     indStepNum = 3;
    indCharge = 4;
    indInitPos = 5:7;
    indFinalPos = 8:10;
%     indTrackLen = 11;
%     indStepLen = 12;
% New TM 
    indInitE = 13;  % New TM (old 14)
    indFinalE = 14; % New TM (old 13) 
    inddE = 15;     % New TM - New
    indEdep = 16;   % New TM  - New Edep that takes over (FinalE - InitE)
    indLayerID = 17;
    
    % N = max(M(:,indTrackID));    %max particle number
    if length(unique(inputMatrix(:,indTrackID)))>1
        error('Multiple particles present in G4 matrix! Please separate before passing to DiffuseTrack4.')
    end
    
    N = size(inputMatrix,1);  %number of steps recorded.
     
    %need to record dE (col 15) for every step written,
    %  and delta-E (col14 - col13) in between consecutive steps of one particle.
   if energyAccountingFlag
        %initalize vectors with enough space for everything (two entries for every row of inputMatrix)
        % intra-step (within step i): index 2*i
        % inter-step (between step i and i+1): index 2*i+1
        depPositions_abs = nan(N*2,3);
        depEnergies = nan(N*2,1);
   else
       %initalize vectors, one for each recorded electron step.
        depPositions_abs = nan(N,3);
        depEnergies = nan(N,1);
   end
    
    if inputMatrix(1,indCharge)==0
        %this particle is a photon. 
        %   photons "deposit" the binding energy of the electron they impact.
        %   this energy will result in an X-ray of <= 1.8 keV (in Si) which will deposit
        %   essentially at the position of the interaction. (very short attenuation length)
        %Any layer skipping, etc. should be handled by parent code.
        if energyAccountingFlag
            depPositions_abs(2:2:2*N,:) = inputMatrix(1:N,indFinalPos); %mm
            %depEnergies(2:2:2*N,1) = inputMatrix(1:N,inddE);            %eV
            depEnergies(2:2:2*N,1) = inputMatrix(1:N,inddE);            %eV
        else
            depPositions_abs(1:1:N,:) = inputMatrix(1:N,indFinalPos); %mm
            depEnergies(1:1:N,1) = inputMatrix(1:N,indEdep);            %eV
        end
            
    elseif inputMatrix(1,indCharge)==-1
        %charged particle, probably an electron.
        if energyAccountingFlag
            %first: intra-step electron deposition. (within a single step)
            %average initial and final positions:
            depPositions_abs(2:2:2*N,:) = (inputMatrix(1:N,indFinalPos) + inputMatrix(1:N,indInitPos))./2; %average initial and final positions
            depEnergies(2:2:2*N,1) = inputMatrix(1:N,indInitE) - inputMatrix(1:N,indFinalE);
            %second: inter-step electron deposition. (between the steps written to file)
            %average final position i and initial position i+1
            depPositions_abs(1:2:2*N-3,:) = (inputMatrix(1:N-1,indFinalPos) + inputMatrix(2:N,indInitPos))./2; 
            depEnergies(1:2:2*N-3,1) = inputMatrix(1:N-1,indFinalE) - inputMatrix(2:N,indInitE);
        else
            %depEnergies(2:2:2*N,1) = inputMatrix(1:N,inddE);
            % Using indEdep instead:
                    % Very simply taking the average position for each
                    % electron step, and using the indEdep for that avg.
                    % position (shouldn't need to do ANYTHING extra for
                    % the first and last positions) 
                    % ...Steps are 1micron steps, start postion shifted by 1/2 micron 
            depEnergies(1:N,1) = inputMatrix(1:N,indEdep);
            depPositions_abs(1:1:N,:) = (inputMatrix(1:N,indFinalPos) + inputMatrix(1:N,indInitPos))./2; 
         
        end
        
    end
    
    %%Unit conversion
    depPositions_abs = depPositions_abs .*1e3;        %mm to um
    depEnergies = depEnergies .*1e-3;     %eV to keV
    
else    %manual mode
    if size(inputMatrix,2) ~=4
        error('Manual inputdata must be 4 columns')
    end
    depPositions_abs(:,1:3) = inputMatrix(:,1:3);
    depEnergies(:,1) = inputMatrix(:,4);
end

% From here on, 'x' references in variable names refer to indLateral(1)
%   'y' references in variable names refer to indLateral(2)

%% Set up image

%Identify and remove points without energy deposition.
%otherwise, empty energy deposition (e.g. from movement of source photon) 
%   might cause the track image to be huge
%this also gets rid of NaN's from the initialization of x and dE 
%   (e.g. from the space for inter-record deposition, where not applicable)
%(these are probably taken care of in Geant4TrackHandling, but just to be safe)
lgNonzeroDeposition = depEnergies>0;
depPositions_abs = depPositions_abs(lgNonzeroDeposition,:);
depEnergies = depEnergies(lgNonzeroDeposition);

xMin = min(depPositions_abs(:, indLateral(1)));
yMin = min(depPositions_abs(:, indLateral(2)));
xMax = max(depPositions_abs(:, indLateral(1)));
yMax = max(depPositions_abs(:, indLateral(2)));

%offsets defined so that (0,0) will be between pixels.
xOffset = floor((xMin - imageEdgeBuffer) / maxPixelSize);   %in pixels
yOffset = floor((yMin - imageEdgeBuffer) / maxPixelSize);   %in pixels

%image size in pixels. 0,0 must be the center of a pixel.
imageSize = [ceil((xMax + imageEdgeBuffer) / maxPixelSize) - xOffset, ...
    ceil((yMax + imageEdgeBuffer) / maxPixelSize) - yOffset];

%adjust positions to be relative to the corner of the image, and the pixel plane. (um still)
depPositions_adj(:,indLateral(1)) = depPositions_abs(:,indLateral(1)) - xOffset*maxPixelSize;
depPositions_adj(:,indLateral(2)) = depPositions_abs(:,indLateral(2)) - yOffset*maxPixelSize;
depPositions_adj(:,indDepth) = depPositions_abs(:,indDepth) - depthCoordinates(1);

%new boundaries on CCD volume
depthCoordinates_adj(1:2) = depthCoordinates(1:2) - depthCoordinates(1);

%initialize fine grid of track image.
imageFine = zeros(imageSize.*maxPixelSize./psftLateralResolution);

%Identify and remove points out of range of the defined CCD volume; 
%   this is not fatal but may indicate bad parameters.
%   Note that this is NOT used for the "contained" flag; that is taken care of in Geant4TrackHandling.
%Both x2 and depthcoordinates are in um now.
lgNotContained = depPositions_adj(:,indDepth) < min(depthCoordinates_adj) | ...
    depPositions_adj(:,indDepth) > max(depthCoordinates_adj);
if any(lgNotContained)
%     warning(['Ignoring ',num2str(sum(lgNotContained)),' points outside of defined CCD geometry; check input parameters.'])
    %ignore points)
    depPositions_adj(lgNotContained,:) = [];
    depEnergies(lgNotContained,:) = [];
end

%% Diffusion.

%lookup psf for each energy deposition point.
%sum into image of resolution 'psftLateralResolution'.
for positionIndex = 1:size(depPositions_adj,1)
    %"abs" accounts for different directions of drift.
    psftIndex = round(abs(depPositions_adj(positionIndex,indDepth)) / psftDepthResolution)+1;
    psfCurrent = depEnergies(positionIndex,1) .* psft{psftIndex};    %(keV)
    
    %round x,y to the nearest fine grid coordinate
    x_psfOffset = round(depPositions_adj(positionIndex,indLateral(1)) / psftLateralResolution);
    y_psfOffset = round(depPositions_adj(positionIndex,indLateral(2)) / psftLateralResolution);
    
    %add psfCurrent onto imageFine
    currentXIndices = (x_psfOffset - psfPatchRadius):(x_psfOffset + psfPatchRadius);
    currentYIndices = (y_psfOffset - psfPatchRadius):(y_psfOffset + psfPatchRadius);
    imageFine(currentXIndices, currentYIndices) = ...
        imageFine(currentXIndices, currentYIndices) + psfCurrent;
end

%% Pixel-size-based data structure.
if length(pixelSize)==1
    %ignore pixel size in the data structure.
    
    %% Image aggregation.
    
    %create CCD image at pixel resolution.
    imageFinal = zeros(imageSize);
    if ~isempty(depEnergies)        %save time if there's nothing here
        for xIndex = 1:imageSize(1)      %x pixel index
            for yIndex = 1:imageSize(2)  %y pixel index
                imageFinal(xIndex,yIndex) = sum(sum(imageFine(...
                    1+(xIndex-1)*pixelSize/psftLateralResolution:xIndex*pixelSize/psftLateralResolution,...
                    1+(yIndex-1)*pixelSize/psftLateralResolution:yIndex*pixelSize/psftLateralResolution)));
            end
        end
    end
    
    % Here is where we would apply the final image orientation, if implemented 
    %   in the future.
    
    %% Outputs
    
    %built cheat structure
    cheat.Edep = sum(depEnergies);     %already in keV
    
    cheat.x = depPositions_abs; %original coordinate system
    cheat.dE = depEnergies;
    
    D.offsets = [xOffset,yOffset];  %in pixels
    D.pixsize = pixelSize;
    D.img = imageFinal;
    D.cheat = cheat;

else
    %need to be careful. And D will be a cell array.
    D = cell(1,length(pixelSize));
    
    for pixelSizeIndex = 1:length(pixelSize)
        currentPixelSize = pixelSize(pixelSizeIndex);
        
        %% Image setup.
        
        %repeat setup, and then simply adjust the fine grid.
        %offsets defined so that (0,0) will be between pixels.
        xOffsetCurrent = floor((xMin - imageEdgeBuffer) / currentPixelSize);   %in pixels
        yOffsetCurrent = floor((yMin - imageEdgeBuffer) / currentPixelSize);   %in pixels
        %it is possible for these offsets to result in image boundaries outside of the 
        %   boundaries for maxPixelSize, if maxPixelSize is not divisible by currentPixelSize.
        %if so, cut into the buffer a bit; this should not be an issue for pixel sizes ~ 50um and less.
        
        %"left" edge (if x is horizontal)
        if xOffsetCurrent*currentPixelSize < xOffset*maxPixelSize
            %shift one pixel: this is sufficient.
            xOffsetCurrent = xOffsetCurrent + 1; %in current pixels
        end
        %"bottom" edge (if y is vertical)
        if yOffsetCurrent*currentPixelSize < yOffset*maxPixelSize
            %shift one pixel: this is sufficient.
            yOffsetCurrent = yOffsetCurrent + 1; %in current pixels
        end
        %"right" edge (if x is horizontal)
        xFarEdge = ceil((xMax + imageEdgeBuffer) / maxPixelSize);
        xFarEdgeCurrent = ceil((xMax + imageEdgeBuffer) / currentPixelSize);
        if xFarEdgeCurrent*currentPixelSize > xFarEdge*maxPixelSize
            %shift one pixel: this is sufficient.
            xFarEdgeCurrent = xFarEdgeCurrent - 1;   %in current pixels
        end
        %"top" edge (if y is vertical)
        yFarEdge = ceil((yMax + imageEdgeBuffer) / maxPixelSize);
        yFarEdgeCurrent = ceil((yMax + imageEdgeBuffer) / currentPixelSize);
        if yFarEdgeCurrent*currentPixelSize > yFarEdge*maxPixelSize
            %shift one pixel: this is sufficient.
            yFarEdgeCurrent = yFarEdgeCurrent - 1;   %in current pixels
        end
        
        %image size in pixels. 0,0 must be the center of a pixel.
        imageSizeCurrent = [xFarEdgeCurrent - xOffsetCurrent, yFarEdgeCurrent - yOffsetCurrent];
        
        %adjust for the different pixel size, to use the same fine grid.
        xFinalOffset = (xOffsetCurrent*currentPixelSize - xOffset*maxPixelSize) / psftLateralResolution;    %in psftLateralResolution
        yFinalOffset = (yOffsetCurrent*currentPixelSize - yOffset*maxPixelSize) / psftLateralResolution;
        %because of boundary checks above, FinalOffsets must be positive.
        %because of pixel size checks above (input adjustment), 
        %   FinalOffsets must be a multiple of psftLateralResolution.
        
        %% Image aggregation.
        
        %create CCD image at pixel resolution.
        imageFinal = zeros(imageSizeCurrent);
        if ~isempty(depEnergies)        %save time if there's nothing here
            for xIndex = 1:imageSizeCurrent(1)      %x pixel index
                for yIndex = 1:imageSizeCurrent(2)  %y pixel index
                    imageFinal(xIndex,yIndex) = sum(sum(imageFine(...
                        xFinalOffset+1+(xIndex-1)*currentPixelSize/psftLateralResolution : xFinalOffset+xIndex*currentPixelSize/psftLateralResolution,...
                        yFinalOffset+1+(yIndex-1)*currentPixelSize/psftLateralResolution : yFinalOffset+yIndex*currentPixelSize/psftLateralResolution)));
                end
            end
        end
        
        % Here is where we would apply the final image orientation, if implemented 
        %   in the future.
        
        %% Outputs
        
        %built cheat structure
        cheat.Edep = sum(depEnergies);     %already in keV
        
        cheat.x = depPositions_abs; %original coordinate system
        cheat.dE = depEnergies;
        
        D{pixelSizeIndex}.offsets = [xOffsetCurrent,yOffsetCurrent];  %in pixels
        D{pixelSizeIndex}.pixsize = pixelSize(pixelSizeIndex);
        D{pixelSizeIndex}.img = imageFinal;
        D{pixelSizeIndex}.cheat = cheat;
    end
end