function out = electronTracksPlotting1(Event, sizeXY, cmaps) 
%Want to Put all of the Electrons of interest Into one 'CCD Image'
%
%
% Coincident Experiment:
% % Number of Pixels per Side:
% sizeXY = [2000 2000];
%
if isempty(sizeXY)
    sizeXY = [2000 2000];
end

% Create Empty CCD 
CCDimg = zeros(sizeXY);

% Define Info of CCDimg
%use real extents of CCD.
pixsize = 10.5;
CCDsize = sizeXY;

            % Assumes origin is at Bottom Left of CCD in lateral dimensions.
            %(x=0,y=0)
            imagebounds = 1e-3*[0, CCDsize(1)*pixsize; ...
                                0, CCDsize(2)*pixsize]; %in mm
                 
            imageboundspix = round(imagebounds *1e3./pixsize)  %in pixels

for i=1:length(Event)
    % Build up Additive Image of CCD
    %img = 
    for j=1:length(Event{1,i}.out.E)
        imgTemp = Event{1,i}.out.T{1,j}.img;
        xyOffsets = [Event{1,i}.out.T{1,j}.x, Event{1,i}.out.T{1,j}.y];
        
        CCDimg = AddToImage(CCDimg, imageboundspix, xyOffsets, imgTemp);
        
    end
end 
%Save Output
out.CCDimg = CCDimg;


% Plot Final  

surf(CCDimg)
%colormap hsv
grid off
axis([ 560 650 1180 1400])

%axis equal

xlabel('X')
ylabel('Y')

%if exist('cmaps')
    colormap(cmaps);
%end


end


%%
function CCDimg = AddToImage(CCDimg,imageboundspix,offsets, imgTemp)
    %truncate at edge of CCD active area
    trunc(1) = max(0, imageboundspix(1,1) - offsets(1));                   %left
    trunc(2) = max(0,-imageboundspix(1,2) + offsets(1) + size(imgTemp,1)); %right
    trunc(3) = max(0, imageboundspix(2,1) - offsets(2));                   %bottom
    trunc(4) = max(0,-imageboundspix(2,2) + offsets(2) + size(imgTemp,2)); %top
    
    x1 = -imageboundspix(1,1) + offsets(1) + (trunc(1)+1:(size(imgTemp,1)-trunc(2)));
    y1 = -imageboundspix(2,1) + offsets(2) + (trunc(3)+1:(size(imgTemp,2)-trunc(4)));
    %add to main image
    CCDimg(x1, y1) = CCDimg(x1,y1) + imgTemp(trunc(1)+1:end-trunc(2), trunc(3)+1:end-trunc(4));
end

