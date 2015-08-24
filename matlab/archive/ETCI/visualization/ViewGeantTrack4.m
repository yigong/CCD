function ViewGeantTrack4(varargin)
%View electron tracks in 3D.
%
% Rewritten for Geant4TrackHandling7 and DiffuseTrack7.
%   x coordinates are absolute, not relative to the image.
%
% ViewGeantTracks(M)            
%    Load geant output matrix M and plot in 3D.
% ViewGeantTracks(M,D)          
%    Load geant output matrix M and DiffuseTrack7 output structure D and plot in 3D.
% ViewGeantTracks(M,D,H,...)
%    Load M and D and HybridTrack algorithm output H. Plot in 3D.
%
%      additional flags:
%
% 'cmap', cmap      specify colormap for track image
% 'renderer', r     specify MATLAB renderer to use
% 'Eplotstyle', s   specify style for plotting energy deposition points
%                   0: plot normal dots
%                   1: vary dot size with energy
%                   2: vary dot color with energy (not working well yet)

%% defaults
diffusionFlag = false;  %use diffusion info
algorithmFlag = false;  %use HybridTrack info
colormapFlag = false;

%% options
%main display options
newFigureFlag = false;  %create new Matlabfigure

initialDirectionTrueFlag = true; %plot true initial electron direction
initialDirectionAlgorithmFlag = true;  %plot initial electron direction as measured by algorithm
algorithmPointsFlag = true;  %plot algorithm points on track image

%Energy-deposited plot styles:
% 0:    plot normal dots
% 1:    vary dot size with energy
% 2:    vary dot color with energy (not working well right now)
energyPlotStyle = 0;

%oldcoordsysflag MUST be set correctly, so VGT knows where the pixel plane is.
%yflipflag is just a visualization preference.
oldCoordinateSystemFlag = true;     %true means depth dimension is 'x' in the geant4 matrix. true for the 13m tracks of ~2010; 
                            %false means depth dimension is 'z' in the geant4 matrix. false for 4/2011 coinc simulation
yFlipFlag = false;          %flip x and y? 'true' -> match TrackGUI images, according to my previous notes..

%graphics customization
rendererflag = true;
renderer = 'zbuffer';

% eplotres = 1;   %condense this many points into one display dat, for eplotstyle>0
energyPlotGain = 7;

zLimitsFlag = true;    %set zlim to 0,650
startMarkerFlag = true; %plot a point at the start of the track

initialDirectionLength = 75; %um in 3D
initialDirectionWidth = 4;   %matlab linewidth units
initialDirectionTrueColor = [0,.8,0];
initialDirectionAlgorithmColor = [0,0.1,.9];

trackLineSpec = '.r';
startLineSpec = ':r';
endLineSpec = ':r';
startMarkerSpec = '*k';
algorithmLineSpec = '.b';
algorithmColor = initialDirectionAlgorithmColor;
algorithmMarkerSize = 10;


pixelSize = 10.5;     %um
buffer = 8.5*pixelSize;    %um

%% Input handling
M = varargin{1};
if nargin==2
    diffusionFlag=true;
    D = varargin{2};
elseif nargin==3
    diffusionFlag=true;
    D = varargin{2};
    if ~isempty(varargin{3}) && isstruct(varargin{3})
        H = varargin{3};
        algorithmFlag = true;
    end
elseif nargin>3
    diffusionFlag = true;
    D = varargin{2};
    if ~isempty(varargin{3}) && isstruct(varargin{3})
        H = varargin{3};
        algorithmFlag = true;
    end
    %look at remaining input arguments
    argind = 4;
    while nargin-argind > 0  %at least two args left to process
        if ischar(varargin{argind}) && strcmpi(varargin{argind},'cmap')
            colormapFlag = true;
            cmap = varargin{argind+1};
        elseif ischar(varargin{argind}) && strcmpi(varargin{argind},'Eplot')
            energyPlotStyle = varargin{argind+1};
        elseif ischar(varargin{argind}) && strcmpi(varargin{argind},'renderer')
            renderer = varargin{argind+1};
        else
            warning('input arguments not recognized...')
        end
        argind = argind + 2;
    end
    
end

%% setup and plotting

if diffusionFlag
    if ~isfield(D,'img')
        return
    end
    %set up diffused track
%     if isfield(D,'cheat') && isfield(D.cheat,'x')
        
%     else
        img = D.img;
        x = D.x2;
        dE = D.dE;
        
        
        %set up image
        minx = min(x(:,1));
        miny = min(x(:,2));
        maxx = max(x(:,1));
        maxy = max(x(:,2));
%     end
end


if newFigureFlag
    if rendererflag
        figure('renderer',renderer); 
        axes('fontsize',14)
    else
        figure;
        axes('fontsize',14)
    end
end
    

%plot the track image
if diffusionFlag
    xi = minx-buffer:pixelSize:minx-buffer+(size(img,1)-1)*pixelSize;
    yi = -miny+buffer:-pixelSize:-miny+buffer-(size(img,2)-1)*pixelSize;
%         x = x(length(x):-1:1);
%         yi = yi(length(yi):-1:1);
    if yFlipFlag
        surf(-yi,xi,img./1e4,'linestyle','none');
    else
        surf(xi,-yi,img'./1e4,'linestyle','none');
    end
    
    if colormapFlag
        colormap(cmap);
    end
    %negative sign flips image to match TrackGUI images
    view(3)
    hold on
end

if oldCoordinateSystemFlag
    %shift coordinate system - for older set of 13m tracks (not april 2011 tracks)
    %this is NOT the same as what yflipflag does!!
    M(:,5:10) = M(:,[6,7,5, 9,10,8]);
end

if yFlipFlag
    %adjust columns of positions in D.x as well as M, so we can forget about it after this
    x(:,1:3) = x(:,[2,1,3]);
    M(:,[5:10]) = M(:,[6,5,7, 9,8,10]);
end

%plot the simulated track points
% x1 = t{i}{1}(:,5).*-1000;
% y1 = t{i}{1}(:,4).*1000;
% z1 = t{i}{1}(:,6).*1000+zp;
% z1 = 650 - z1;

if energyPlotStyle==0
    plot3(x(:,1),x(:,2),x(:,3),trackLineSpec)
elseif energyPlotStyle==1
    for i=1:size(x,1)
        plot3(x(i,1),x(i,2),x(i,3),trackLineSpec,'markersize',dE(i)*energyPlotGain);
    end
elseif energyPlotStyle==2
    dEmax = max(dE);
    minGrayLevel = 0.8; %this is "zero"
    for i=1:size(x,1)
        %'color' argument overrides linespec
        plot3(x(i,1),x(i,2),x(i,3),trackLineSpec,'color',ones(1,3)*(minGrayLevel-(minGrayLevel)*dE(i)/dEmax),'markersize',3);
    end
else
    
end
hold on

if diffusionFlag
    %plot nice dashed lines at start and end
    primaryind = M(find(M(:,2)==1,1),1);
    x0 = M(find(M(:,1)==primaryind,1),5:7) *1e3;    %mm to um
    xf = M(find(M(:,1)==primaryind,1,'last'),8:10) *1e3;
    %need to translate to the coordinate system of the plot..
    %use the knowledge that x0 and x2(1,:) should be the same physical location
    dx(1:2) = x(1,1:2) - x0(1:2);
    x0(1:2) = x0(1:2) + dx(1:2);
    xf(1:2) = xf(1:2) + dx(1:2);
    x0(3) = abs(x0(3));
    xf(3) = abs(xf(3));
    plot3([x0(1),x0(1)],[x0(2),x0(2)],[0,x0(3)],startLineSpec,'linewidth',2)
    plot3([xf(1),xf(1)],[xf(2),xf(2)],[0,xf(3)],endLineSpec,'linewidth',2)
    
    if startMarkerFlag
        plot3(x0(1),x0(2),x0(3),startMarkerSpec)
    end
    
    if initialDirectionTrueFlag
        %plot initial direction of electron
        sina = sin(D.cheat.a*pi/180);
        cosa = cos(D.cheat.a*pi/180);
        sinb = sin(D.cheat.b*pi/180);
        cosb = cos(D.cheat.b*pi/180);
        initdirtruex = [x0(1),x0(1) + initialDirectionLength*cosa*cosb];
        initdirtruey = [x0(2),x0(2) + initialDirectionLength*sina*cosb];
        initdirtruez = [x0(3),x0(3) + initialDirectionLength*sinb];
        if yFlipFlag
            plot3(initdirtruey,initdirtruex,initdirtruez,'-','linewidth',initialDirectionWidth,'color',initialDirectionTrueColor);
        else
            plot3(initdirtruex,initdirtruey,initdirtruez,'-','linewidth',initialDirectionWidth,'color',initialDirectionTrueColor);
        end
        
    end
end

if algorithmFlag
    if algorithmPointsFlag && isfield(H,'x')
        %plot algorithm points
        if yFlipFlag
            xpoints = (H.y - 8) * pixelSize;
            ypoints = (H.x - 8) * pixelSize;
        else
            xpoints = (H.x - 8) * pixelSize;
            ypoints = (H.y - 8) * pixelSize;
        end
        zpoints = 1e-2 * ones(size(xpoints));
        plot3(xpoints,ypoints,zpoints,algorithmLineSpec,'color',algorithmColor,'markersize',algorithmMarkerSize)
    end
    
    if initialDirectionAlgorithmFlag && isfield(H, 'alpha')
        %plot initial direction of electron, measured by algorithm
        sina = sin(H.alpha*pi/180);
        cosa = cos(H.alpha*pi/180);
        sinb = sin(H.beta*pi/180);
        cosb = cos(H.beta*pi/180);
        %plot both plus and minus beta
        initdiralgx = [x0(1) + initialDirectionLength*cosa*cosb, x0(1), x0(1) + initialDirectionLength*cosa*cosb];
        initdiralgy = [x0(2) + initialDirectionLength*sina*cosb, x0(2), x0(2) + initialDirectionLength*sina*cosb];
        initdiralgz = [x0(3) + initialDirectionLength*sinb, x0(3), x0(3) - initialDirectionLength*sinb];
        if yFlipFlag
            plot3(initdiralgy,initdiralgx,initdiralgz,'-','linewidth',initialDirectionWidth,'color',initialDirectionAlgorithmColor);
        else
            plot3(initdiralgx,initdiralgy,initdiralgz,'-','linewidth',initialDirectionWidth,'color',initialDirectionAlgorithmColor);
        end
    end
end

axis equal
if zLimitsFlag
    zlim([0,650]);
end

grid off
hold off
