function [D] = DiffuseTrack5(varargin)
%function [D] = DiffuseTrack4(psft,mode,inputdata)
%function [D] = DiffuseTrack4(psft,mode,inputdata,'thickness',th)
%function [D] = DiffuseTrack4(psft,mode,inputdata,'plotflag')
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
%   thickness: thickness of CCD, in um. (default 650)
%   plotflag: sets plotflag to true. (default false)
%   zdepth: sets zdepthflag to true, so columns 7 & 10 are the depth dimension. (default false)
%   negdepth: sets znegflag to true, so the back face is at -th/1000 instead of +th/1000.
%   depthinvert: sets zinvertflag to true, so the back face is at 0 and the 
%       pixel plane is at +th/1000 (or -th/1000 if negdepth)
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
%     cheat.x0:    x coordinate of start of trajectory (pixels)
%     cheat.y0:    y coordinate of start of trajectory (pixels)
%     cheat.x:       list of 3D coordinates that have been diffused. z=depth.
%     cheat.xp:      list of 2D coordinates, in units of pixels from corner of image.
%     cheat.dE:      list of energies diffused. corresponds to D.x.

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
    M = varargin{3};
else
    error('inputdata should be a numeric array')
end

%defaults
th = 650;
plotflag = false;
zdepthflag = false;
znegflag = false;
zinvertflag = false;
energyaccountingflag = false;
everystepflag = false;

%non-defaults
i=4;
while true
    if i > nargin
        break
    end
    if strcmpi(varargin{i},'thickness') || strcmpi(varargin{i},'th')
        if length(varargin)>i
            if isnumeric(varargin{i+1})
                th = varargin{i+1};
                i=i+2;
            else
                error('thickness should be numeric')
            end
        else
            warning('optional argument "thickness" not defined')
            i=i+1;
        end
    elseif strcmpi(varargin{i},'plotflag') || strcmpi(varargin{i},'plot')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==0 || varargin{i+1}==1)))
            plotflag = logical(varargin{i+1});
            i=i+2;
        else
            plotflag = true;
            i=i+1;
        end
    elseif strcmpi(varargin{i},'zdepthflag') || strcmpi(varargin{i},'zdepth') || strcmpi(varargin{i},'zdep')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==0 || varargin{i+1}==1)))
            zdepthflag = logical(varargin{i+1});
            i=i+2;
        else
            zdepthflag = true;
            i=i+1;
        end
    elseif strcmpi(varargin{i},'znegflag') || strcmpi(varargin{i},'negdepth') || strcmpi(varargin{i},'znegative')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==0 || varargin{i+1}==1)))
            znegflag = logical(varargin{i+1});
            i=i+2;
        else
            znegflag = true;
            i=i+1;
        end
    elseif strcmpi(varargin{i},'zinvertflag') || strcmpi(varargin{i},'zinv') || strcmpi(varargin{i},'depthinvert')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==0 || varargin{i+1}==1)))
            zinvertflag = logical(varargin{i+1});
            i=i+2;
        else
            zinvertflag = true;
            i=i+1;
        end
    elseif strcmpi(varargin{i},'energyaccountingflag') || strcmpi(varargin{i},'energyaccounting') || strcmpi(varargin{i},'manualenergycorrect')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==0 || varargin{i+1}==1)))
            energyaccountingflag = logical(varargin{i+1});
            i=i+2;
        else
            energyaccountingflag = true;
            i=i+1;
        end
    end
end


%
% New in DiffuseTrack3:
%  geant4 must write any particle step that produces a secondary.
%    This means that secondary electrons need not be matched to their parent,
%    because there is no ambiguity in deposited energy.
%  multiple electrons (from a single photon parent) are all tracked, and output
%    as cells.

pixsize = 10.5;     %um
buffer = 8.5*pixsize;    %um            %point (0,0) is in the center of a pixel

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
    if length(unique(M(:,1)))>1
        error('Multiple particles present in G4 matrix! Please separate before passing to DiffuseTrack4.')
    end
    
    N = size(M,1);  %number of steps recorded.
    
    Etot = M(1,indInitE) .* 1e-3; %first E_initial, in keV
    %%Step info.
    
    %need to record dE (col 15) for every step written,
    %  and delta-E (col14 - col13) in between consecutive steps of one particle.
    
    %define variables: 
    %   x (position of energy deposition, (mm,mm,mm))
    %   dE (amount of energy deposited, eV)
    %   ss (step size.. do I care?)
    
    %initalize vectors with enough space for everything (two entries for every row of M)
    % intra-step (within step i): index 2*i
    % inter-step (between step i and i+1): index 2*i+1
    x = nan(size(M,1)*2,3);
    dE = nan(size(M,1)*2,1);
    % ss = nan(size(M,1)*2,1);
    
    %first: photons "deposit" the binding energy of the electron they impact.
    %       this energy will result in an X-ray of <= 1.8 keV (in Si) which will deposit
    %       essentially at the position of the interaction. (very short attenuation length)
    if M(1,indCharge)==0
        %this particle is a photon. Any layer skipping, etc. should be handled by parent code.
        x(2:2:2*N,:) = M(1:N,indFinalPos);
        dE(2:2:2*N,1) = M(1:N,inddE);
    end
    
    if M(1,indCharge)==-1
        %charged particle, probably an electron.
        
        %first: intra-step electron deposition.
        x(2:2:2*N,:) = (M(1:N,indFinalPos) + M(1:N,indInitPos))./2; %average initial and final positions
        if energyaccountingflag
            dE(2:2:2*N,1) = M(1:N,indInitE) - M(1:N,indFinalE);
        else
            dE(2:2:2*N,1) = M(1:N,inddE);
        end
        
        %second: inter-step electron deposition.
        x(1:2:2*N-3,:) = (M(1:N-1,indFinalPos) + M(2:N,indInitPos))./2; %average final position i and initial position i+1
        dE(1:2:2*N-3,1) = M(1:N-1,indFinalE) - M(2:N,indInitE);
    end
    
    %%Unit conversion, coordinate transform.
    
    x = x .*1e3;        %mm to um
    dE = dE .*1e-3;     %eV to keV
    
    if ~zdepthflag
        %  x(:,1) is the depth dimension. ('x' => 'z')
        %  x(:,2) is the first pixel dimension. ('y' => 'x')
        %  x(:,3) is the second pixel dimension. ('z' => 'y')
        x2(:,3) = x(:,1);
        x2(:,1) = x(:,2);
        x2(:,2) = x(:,3);
        x = x2; clear x2
    end
    
    %get CCD thickness onto interval [0,+650] according to user input
    if znegflag
        x(:,3) = -x(:,3);
    end
    
    if zinvertflag
        x(:,3) = 650 - x(:,3);
    end
    
else    %manual mode
    if size(M,2) ~=4
        error('Manual inputdata must be 4 columns')
    end
    x(:,1:3) = M(:,1:3);
    dE(:,1) = M(:,4);
    Etot = sum(dE);
end

x = x(dE>0,:);
%without dE>0, empty energy deposition (e.g. from movement of source photon) might cause the track image to be huge
%this also gets rid of NaN's from the initialization of x and dE (e.g. from the space for inter-record deposition, where not applicable)
dE = dE(dE>0);

%% Set up image

% zres = 0.5; %um resolution of PSFT
zres = th/(length(psft)-1);

int=0.5; %um
%This is the resolution of the fine grid. Energy deposition coordinates
%will be rounded off to a grid position.

%int must match the xres of PSFT

minx=min(x(:,1));
miny=min(x(:,2));
maxx=max(x(:,1));
maxy=max(x(:,2));

xoffset = floor((minx - buffer)/pixsize);   %in pixels
yoffset = floor((miny - buffer)/pixsize);   %in pixels
%imsize in pixels. 0,0 must be the center of a pixel.
imsize = [ceil((maxx+buffer)/pixsize) - xoffset, ...
    ceil((maxy+buffer)/pixsize) - yoffset];

x2(:,1) = x(:,1) - xoffset*pixsize;  %microns from image edge.
x2(:,2) = x(:,2) - yoffset*pixsize;  %microns from image edge.
x2(:,3) = x(:,3);   %depth
%x2 is defined such that (0,0) is a pixel edge/corner.

%initialize fine grid of track image.
fgrid=zeros(imsize.*pixsize./int);

gs = size(psft{1},1);    %dimensions of psf.
gs2 = (gs-1)/2;

%% Diffusion.

%lookup psf for each energy deposition point.
%sum into image of resolution 'int'.
for i=1:size(x2,1)
    zi = round(x2(i,3)/zres)+1;     %index of psft to use
    Espread = dE(i,1).*psft{zi};    %(already keV units)
    xi = round(x2(i,1)/int);        %round x,y to the nearest fine grid coordinate
    yi = round(x2(i,2)/int);
    fgrid(xi-gs2:xi+gs2,yi-gs2:yi+gs2) = fgrid(xi-gs2:xi+gs2,yi-gs2:yi+gs2) + Espread;
end

%% Image aggregation.

%create CCD image at pixel resolution.
img=zeros(imsize);
for p=1:imsize(1)
    for r=1:imsize(2)
%         SGrid(p,r)=sum(sum(Grid(51+pxsize*(p-1)/int:50+pxsize*(p)/int,51+pysize*(r-1)/int:50+pysize*(r)/int)));  
        img(p,r) = sum(sum(fgrid(1+(p-1)*pixsize/int:p*pixsize/int,1+(r-1)*pixsize/int:r*pixsize/int)));
    end
end

%% Outputs.

cheat.Edep = sum(dE);     %already in keV
cheat.Etot = Etot; %first E_initial, in keV
cheat.x0 = x2(1,1:3);  %in um, image coordinates

cheat.x = x2;
cheat.xp = x2(:,1:2)./pixsize;
cheat.dE = dE;

D.offsets = [xoffset,yoffset];  %in pixels
D.pixsize = pixsize;

D.img = img;
D.cheat = cheat;

if plotflag
    surf(img,'linestyle','none')
    view(2)
    colormap('hot')
    axis equal
end

