function h = SurfElectronTrack(varargin)
% function h = SurfElectronTrack(img);
% function h = SurfElectronTrack(img,'title',plot_title);
% function h = SurfElectronTrack(img,'cmap',cmap);
% function h = SurfElectronTrack(img,'FixAxesMargins');
% function h = SurfElectronTrack(img,'h',axes_handle);
% function h = SurfElectronTrack(img,'clim',[color_min, color_max]);
% function h = SurfElectronTrack(img,'clim','auto');
% function h = SurfElectronTrack(img,'x',[x',y']);
% function h = SurfElectronTrack(img,'offsets',[xoffset,yoffset]);
%
% Make a 2D surface plot of an electron track.
%
% X label and Y label are set already.
% img must be supplied as an array of numerical values.
% Title is blank unless supplied as an argument. (string)
% Colormap is 'jet' unless supplied as an argument. (Should use cmaphotlog.)
% h plots the track onto existing axes with handle h.
% clim sets the range of the color axis to [color_min, color_max].

if nargin<1 || ~isnumeric(varargin{1})
    error('Need numeric img input')
end
img = varargin{1};

%defaults
plot_title = [];
xPoints = [];
cmap = 'jet';
xLabel = 'Pixels (10.5 \mum)';
yLabel = xLabel;
FixAxesMarginsFlag = false;
NewAxesFlag = true;     %a bit redundant but ok whatevers
h = [];
pixelSize = 10.5;
ColorLimits = [];

%input handling
n=2;
while n<nargin
    if strcmpi(varargin{n},'title')
        plot_title = varargin{n+1};
        n=n+2;
    elseif strcmpi(varargin{n},'cmap') || strcmpi(varargin{n},'colormap')
        cmap = varargin{n+1};
        n=n+2;
    elseif strcmpi(varargin{n},'FixAxesMargin') || strcmpi(varargin{n},'FixAxesMargins')
        FixAxesMarginsFlag = true;
        n=n+1;
    elseif strcmpi(varargin{n},'h') || strcmpi(varargin{n},'handle')
        NewAxesFlag = false;
        h = varargin{n+1};
        n=n+2;
    elseif strcmpi(varargin{n},'clim')
        if isnumeric(varargin{n+1})
            if length(varargin{n+1}) ~= 2
                error('Color limits should be a 2-element vector.')
            end
            ColorLimits = varargin{n+1};
        else
            if ischar(varargin{n+1}) && strcmpi(varargin{n+1},'auto')
                ColorLimits = [];
            else
                error('Color limits should be 2-element vector, or ''auto''.')
            end
        end
        n=n+2;
    elseif strcmpi(varargin{n},'x')
        if isnumeric(varargin{n+1})
            xPoints = varargin{n+1};
            n=n+2;
        else
            error('x should be a numeric array')
        end
    elseif strcmpi(varargin{n},'offsets') && length(varargin{n+1})==2
        offsets = varargin{n+1};
        n=n+2;
    elseif strcmpi(varargin{n},'pixelsize') || strcmpi(varargin{n},'pixsize')
        pixelSize = varargin{n+1};
        
        %this stuff shouldn't be here
        xLabel = ['Pixels (',num2str(pixelSize),' \mum)'];
        yLabel = xLabel;
        n=n+2;
    else
        error('I don''t understand your argument ID.')
    end
end

% ColorLimits = [0,25*pixelSize^2/10.5^2];    %scale to pixel size

%add buffer on far edge
img = [img, zeros(size(img,1),1); zeros(1,size(img,2)+1)];

if NewAxesFlag
    %generate plot
    figure('color','w');
    % hf = figure('color','none');
    h = axes('fontsize',18);
else
    %bring figure/axes to focus
    hf = get(h,'Parent');
    figure(hf);
    set(hf,'CurrentAxes',h);
    hold off;
end
pcolor(img); 
shading('flat');
view(2); 
axis equal tight; 
colormap(cmap); 
colorbar; 
if ~isempty(ColorLimits)
    caxis(ColorLimits);
end
if ~isempty(plot_title>0)
    title(plot_title,'fontsize',18); 
end
if ~isempty(xPoints)
    hold on;
    %invert x and y for pcolor/surf
    plot3(xPoints(:,2)/pixelSize - offsets(2) + 1, xPoints(:,1)/pixelSize - offsets(1) + 1, ones(size(xPoints,1),1), '.c');
end
set(h,'TickDir','out');
xlabel(xLabel,'fontsize',18);
ylabel(yLabel,'fontsize',18);
drawnow;
if FixAxesMarginsFlag
    FixAxesMargins(h);
end
