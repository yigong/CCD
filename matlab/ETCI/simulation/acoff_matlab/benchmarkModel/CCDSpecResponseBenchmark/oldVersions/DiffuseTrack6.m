function [D] = DiffuseTrack6(varargin)
%function [D] = DiffuseTrack6(psft, mode, inputdata)
%function [D] = DiffuseTrack6(psft, mode, inputdata, argname, argvalue, ...)
%function [D] = DiffuseTrack6(psft, mode, inputdata, argflagname, ...)
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
%   CCDx: defines how the image is oriented in data space.
%       
%   pixelsize: pixel pitch, in um. (default 10.5)
%   energyaccounting: for use with Geant4TrackHandling, manual energy correction mode. 
%       uses InitE - FinalE, instead of dE, for intra-step deposition.
%
%output variables
%
%D is a cell array, with one cell for each primary electron in this event sequence.
%
%D{i}.img:     data matrix of ccd plane pixels.    (keV)
%D{i}.offsets: absolute position of pixel (i,j) is (offsets(1)+i,offsets(2)+j). origin matches that of input data. 
%D{i}.cheat:
%     cheat.Etot:    actual electron energy        (keV)
%     cheat.Edep:    energy deposited in detector  (keV)
%     cheat.z0:      zp (starting z)               (um)
%     cheat.x0:      x coordinate of start of trajectory (pixels)
%     cheat.y0:      y coordinate of start of trajectory (pixels)
%     cheat.x:       list of 3D coordinates that have been diffused. z=depth.
%     cheat.xp:      list of 2D coordinates, in units of pixels from corner of image.
%     cheat.dE:      list of energies diffused. corresponds to D.x.
%
% future work: 
% change format of psft to include thickness or depth resolution.
%   currently it will be hardcoded to 650um thickness, so that we can apply it to 
%   a "thinner" but actually non-depleted CCD.
% make the fine grid resolution an optional argument (also depends on psft)
% review outputs
% should cheat alpha, beta change if the track does not start inside the CCD geometry?
%   we shouldn't encourage or allow this case... should be dealt with in Geant4TrackHandling


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
    if strcmpi(varargin{i},'pixelsize') || strcmpi(varargin{i},'pixelpitch') || strcmpi(varargin{i},'pixel')
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
    elseif strcmpi(varargin{i},'driftdim') || strcmpi(varargin{i},'drift')
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
    end
end

%% Fix up parameters.

%Convert depthcoordinates from mm to um
depthCoordinates = depthCoordinates*1e3;

psfPatchWidth = size(psft{1},1);    %dimensions of psf, in units of the fine grid
psfPatchRadius = (psfPatchWidth-1)/2;

%Check for geometry errors
if (abs(depthCoordinates(2) - depthCoordinates(1)) - psftThickness) > 0    %allow rounding error?
    error('CCD geometry should not be thicker than CCD thickness described by diffusion table.')
end

%space around the energy deposition points, to accomodate plenty of diffusion.
%   this buffer is enough for any pixel size, don't mind the hardcoded 10.5.
imageEdgeBuffer = 8.5*10.5;    %um            %point (0,0) is in the center of a pixel

%% Set up energy deposition points.

if g4mode
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
    
    % N = max(M(:,indTrackID));    %max particle number
    if length(unique(inputMatrix(:,indTrackID)))>1
        error('Multiple particles present in G4 matrix! Please separate before passing to DiffuseTrack4.')
    end
    
    N = size(inputMatrix,1);  %number of steps recorded.
    
    totalEnergy = inputMatrix(1,indInitE) .* 1e-3; %first E_initial, in keV (from eV)
     
    %need to record dE (col 15) for every step written,
    %  and delta-E (col14 - col13) in between consecutive steps of one particle.
    
    %initalize vectors with enough space for everything (two entries for every row of inputMatrix)
    % intra-step (within step i): index 2*i
    % inter-step (between step i and i+1): index 2*i+1
    depPositions_abs = nan(N*2,3);
    depEnergies = nan(N*2,1);
    
    if inputMatrix(1,indCharge)==0
        %this particle is a photon. 
        %   photons "deposit" the binding energy of the electron they impact.
        %   this energy will result in an X-ray of <= 1.8 keV (in Si) which will deposit
        %   essentially at the position of the interaction. (very short attenuation length)
        %Any layer skipping, etc. should be handled by parent code.
        depPositions_abs(2:2:2*N,:) = inputMatrix(1:N,indFinalPos); %mm
        depEnergies(2:2:2*N,1) = inputMatrix(1:N,inddE);            %eV
    elseif inputMatrix(1,indCharge)==-1
        %charged particle, probably an electron.
        
        %first: intra-step electron deposition. (within a single step)
        depPositions_abs(2:2:2*N,:) = (inputMatrix(1:N,indFinalPos) + inputMatrix(1:N,indInitPos))./2; %average initial and final positions
        if energyAccountingFlag
            depEnergies(2:2:2*N,1) = inputMatrix(1:N,indInitE) - inputMatrix(1:N,indFinalE);
        else
            depEnergies(2:2:2*N,1) = inputMatrix(1:N,inddE);
        end
        
        %second: inter-step electron deposition. (between the steps written to file)
        %average final position i and initial position i+1
        depPositions_abs(1:2:2*N-3,:) = (inputMatrix(1:N-1,indFinalPos) + inputMatrix(2:N,indInitPos))./2; 
        depEnergies(1:2:2*N-3,1) = inputMatrix(1:N-1,indFinalE) - inputMatrix(2:N,indInitE);
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
    totalEnergy = sum(depEnergies);
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
lgNoDeposition = depEnergies>0;
depPositions_abs = depPositions_abs(lgNoDeposition,:);
depEnergies = depEnergies(lgNoDeposition);

xMin = min(depPositions_abs(:, indLateral(1)));
yMin = min(depPositions_abs(:, indLateral(2)));
xMax = max(depPositions_abs(:, indLateral(1)));
yMax = max(depPositions_abs(:, indLateral(2)));

%offsets defined so that (0,0) will be between pixels.
xOffset = floor((xMin - imageEdgeBuffer) / pixelSize);   %in pixels
yOffset = floor((yMin - imageEdgeBuffer) / pixelSize);   %in pixels

%image size in pixels. 0,0 must be the center of a pixel.
imageSize = [ceil((xMax + imageEdgeBuffer) / pixelSize) - xOffset, ...
    ceil((yMax + imageEdgeBuffer) / pixelSize) - yOffset];

%adjust positions to be relative to the corner of the image, and the pixel plane. (um still)
depPositions_adj(:,indLateral(1)) = depPositions_abs(:,indLateral(1)) - xOffset*pixelSize;
depPositions_adj(:,indLateral(2)) = depPositions_abs(:,indLateral(2)) - yOffset*pixelSize;
depPositions_adj(:,indDepth) = depPositions_abs(:,indDepth) - depthCoordinates(1);

%initialize fine grid of track image.
imageFine = zeros(imageSize.*pixelSize./psftLateralResolution);

%Identify and remove points out of range of the defined CCD volume; 
%   this is not fatal but may indicate bad parameters.
%   Note that this is NOT used for the "contained" flag; that is taken care of in Geant4TrackHandling.
%Both x2 and depthcoordinates are in um now.
lgNotContained = depPositions_adj(:,indDepth) < min(depthCoordinates) | ...
    depPositions_adj(:,indDepth) > max(depthCoordinates);
if any(lgNotContained)
    warning(['Ignoring ',num2str(sum(lgNotContained)),' points outside of defined CCD geometry; check input parameters.'])
    %ignore points)
    depPositions_adj(lgNotContained,:) = [];
    depEnergies(lgNotContained,:) = [];
end

%% Diffusion.

%lookup psf for each energy deposition point.
%sum into image of resolution 'psftLateralResolution'.
for i = 1:size(depPositions_adj,1)
    %"abs" accounts for different directions of drift.
    psftIndex = round(abs(depPositions_adj(i,indDepth)) / psftDepthResolution)+1;
    psfCurrent = depEnergies(i,1) .* psft{psftIndex};    %(keV)
    
    %round x,y to the nearest fine grid coordinate
    x_psfOffset = round(depPositions_adj(i,indLateral(1)) / psftLateralResolution);
    y_psfOffset = round(depPositions_adj(i,indLateral(2)) / psftLateralResolution);
    
    %add psfCurrent onto imageFine
    currentXIndices = (x_psfOffset - psfPatchRadius):(x_psfOffset + psfPatchRadius);
    currentYIndices = (y_psfOffset - psfPatchRadius):(y_psfOffset + psfPatchRadius);
    imageFine(currentXIndices, currentYIndices) = ...
        imageFine(currentXIndices, currentYIndices) + psfCurrent;
end

%% Image aggregation.

%create CCD image at pixel resolution.
imageFinal = zeros(imageSize);
for p = 1:imageSize(1)      %x pixel index
    for r = 1:imageSize(2)  %y pixel index
        imageFinal(p,r) = sum(sum(imageFine(...
            1+(p-1)*pixelSize/psftLateralResolution:p*pixelSize/psftLateralResolution,...
            1+(r-1)*pixelSize/psftLateralResolution:r*pixelSize/psftLateralResolution)));
    end
end

% Here is where we would apply the final image orientation, if implemented 
%   in the future.

%% Outputs

%built cheat structure
cheat.Edep = sum(depEnergies);     %already in keV
cheat.Etot = totalEnergy; %first E_initial, in keV
cheat.x0 = depPositions_adj(1,[indLateral(1),indLateral(2),indDepth]);  %in um, image coordinates

cheat.x = depPositions_adj; %note that these positions remain in the original coordinate system
cheat.xp = depPositions_adj(:,1:2)./pixelSize;  %lateral (2D) positions of points, in units of pixels
cheat.dE = depEnergies;

if length(depEnergies)<2
    warning('Only one point of energy deposition: no direction calculated')
else
    dx = cheat.x(2,:) - cheat.x(1,:);
    dl = sqrt(sum(dx(1:3).^2));
    
    %alpha (degrees) measured from indLateral(1) toward indLateral(2)
    cheat.a = atan2(dx(indLateral(2)),dx(indLateral(1))) *180/pi;
    %beta positive is away from the pixel plane
    if depthCoordinates(1) < depthCoordinates(2)
        cheat.b = asin(dx(indDepth)/dl) *180/pi;
    else
        cheat.b = asin(-dx(indDepth)/dl) *180/pi;
    end
    
end


D.offsets = [xOffset,yOffset];  %in pixels
D.pixsize = pixelSize;
D.img = imageFinal;
D.cheat = cheat;