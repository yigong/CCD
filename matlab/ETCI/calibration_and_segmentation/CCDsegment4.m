function [T,E,maskImage] = CCDsegment4(varargin)
% function [T,E] = CCDsegment4(img,th)
% function [T,E] = CCDsegment4(img,th,'Mode',modestr)
% function [T,E] = CCDsegment4(img,th,'PixelThreshold',NumPixelsThreshold)
% function [T,E] = CCDsegment4(img,th,'NeighborLayers',NeighborLayers)
% function [T,E] = CCDsegment4(img,th,'EdgeVeto',edgeVeto)
% function [T,E] = CCDsegment4(img,th,'Subarray',[w,h,colc,rowc])
% function [T,E] = CCDsegment4(img,th,'pixelSize',pixelSize)
% function [T,E] = CCDsegment4(img,th,'UseBWLabel')
% function [T,E,mask] = CCDsegment4(...,'MaskFlag',maskFlag)
%
%Requires the Image Processing toolbox.
% 12/4/12, CCDsegment4: now adjusting neighbor behavior for small pixels.
%
%Inputs:
% img: input image array; should come from CCDcal_median or similar.
%   The black level should be normalized to 0. In 'precal' mode, the image
%   should be in detector units (gain = 1); in 'postcal' mode, the image should
%   be in keV.
% th: energy threshold for segmentation. Units defined by mode.
%     Precal mode: units are [sigma of black level noise]; try 6 or so.
%     Postcal mode: units are keV. Try 0.1 for normal images.
%       **In April 2010 I showed 0.1 keV on my ARI poster, to achieve the 0.571 keV Am-241 resolution.
%         This is 5.3 sigma of standard noise.
%       In October 2012 on my IEEE poster, I reached a limit of around 3.3 sigma noise for smaller pixels
%         and higher noise levels. But this was very rough.
%
% modestr: either 'precal' or 'postcal'
%   'Precal' indicates that the gain of the device was unknown, so the images
%   are in detector units, and are not necessarily the same from quadrant to quadrant.
%   'Postcal' (default) means that the image is calibrated in keV.
% pixth: number indicating minimum number of pixels allowed in a segmented track.
%   If a track is found with n_pixels < pixth, it will not be stored.
% NeighborLayers: The number of layers of neighboring pixels to add to the track.
%   Added with bwmorph 'dilate' operation.
%   One layer of neighbors has been found to improve energy resolution
%   in the Am-241 photopeak for experimental tracks (10.5um, 19eV sigma noise).
%*NeighborFlag operation is still allowed, for backwards compatibility. (sets NeighborLayers=1)
% edgeveto: logical; true means that tracks which touch the edge of the image will be vetoed.
% subarray: [w,h,colc,rowc] parameters to compensate for subarray readout of the CCD. 
%   Use if the Subarray parameters were set in the acquisition software.
%   The parameters are used in CCDsegment3 to label x and y coordinates correctly
%   (in the coordinate system of the full device).
%     w = box width, h = box height, colc = box column center, rowc = box rows center
%     The bias width is assumed to have been set to zero
% MaskFlag: also output a logical array showing mask, including pixels above 
%   the low threshold and their neighbors (if NeighborFlag).
% maxSegments: if more than maxSegments different segments are found,
%   abort. (for smaller pixels and noisier images, when a large image might give 5000 segments)
% pixelSize: in um. Determines the number of layers of neighboring pixels to use. max(1,round(pixelSize/10.5))
%   If NeighborLayers and NumPixelsThreshold are manually defined, then pixelSize is unnecessary.
% UseBWLabel: use CCDsegment3's bwlabel call instead of the new, better bwconncomp. 
% 
% Default options:
%   Mode = 'Postcal'
%   PixelThreshold = 4
%   NeighborLayers = 1
%   EdgeVeto = true
%   MaskFlag = false
%   maxSegments = 500
%
%Outputs:
% T: cell array of structures.
%   T{i}.img: image of track i.
%   T{i}.E: total track energy = sum(T{i}.img(:))
%   T{i}.x: The x offset of the track, in pixels.
%       Here, x is the first array index; +x is up in physical space.
%       The first row in the track image, T{i}.img(1,:), 
%       is the (x+1)th row in the full device image, img(x,:).
%   T{i}.y: The y offset of the track, in pixels.
%       Here, y is the second array index; +y is ... (need to double check which direction) in physical space.
%       The first column in the track image, T{i}.img(:,1),
%       is the (y+1)th column in the full device image, img(y,:).
%   T{i}.edgeflag: marks if the track was neighboring the 
%       edge regions blacked out by CCDcal.
% E: list of energies. Either by quadrant (precal mode) or for the entire image.
% mask: binary mask of track pixels & neighbors.

%input handling
[imageOriginal, threshold, opts] = CCDsegment4_InputHandling(nargin, varargin(:));
    
if opts.maskFlag && opts.edgeVeto
    error('Maskflag and Edgeveto should not be used together...')
    %why? no good reason. These flags are old and outdated I think. #futurework
end


%% 
maskImage = false(size(imageOriginal));
%measure dimensions from the CCD image given
[dim1,dim2] = size(imageOriginal);

if opts.quadrantStraddleVeto
    %Segment each quadrant separately.
    if opts.preCalMode
        %list energy values separately by quadrant, not all in one list.
        %   because each quadrant will have a slightly different gain value.
        E = cell(1,4);
        T = cell(1,4);  %for each quadrant
        
        for quadrantIndex = 1:4     %each quadrant
            %get quadrant image and indices
            [xOffset, yOffset, imageQuadrantOriginal] = ...
                CCDsegment4_GetQuadrant(dim1, dim2, imageOriginal, quadrantIndex);
            
            %smooth, if desired
            if opts.useSmoothing
                imageQuadrant = CCDsegment4_Smooth(opts, imageQuadrantOriginal);
            else
                imageQuadrant = imageQuadrantOriginal;
            end
            
            %determine threshold based on noise level
            threshold = CCDsegment4_GetPreCalThreshold(imageQuadrant,threshold);
            
            if opts.useBWLabel
                [T{quadrantIndex}, E{quadrantIndex}, ind, maskImage] = ...
                    CCDsegment4_bwlabel(imageQuadrant, imageQuadrantOriginal, xOffset, yOffset, threshold, ...
                    opts, T{quadrantIndex}, E{quadrantIndex}, 1);
            else
                [T{quadrantIndex}, E{quadrantIndex}, ind, maskImage] = ...
                    CCDsegment4_bwconncomp(imageQuadrant, imageQuadrantOriginal, xOffset, yOffset, threshold, ...
                    opts, T{quadrantIndex}, E{quadrantIndex}, 1);
            end
            
            %cut off any extra cells or NaN's from T{quadrantIndex} and E{quadrantIndex}.
            T{quadrantIndex} = T{quadrantIndex}(1:ind-1);
            E{quadrantIndex} = E{quadrantIndex}(1:ind-1);
        end        
    else    % [CCDsegment4: currently quadrantStraddleVeto==PreCalMode so this doesn't ever happen]
        %energy values from all quadrants will all be in one list
        E = [];
        T = cell(0,1);
        %keep an index of where we are in T and E
        ind = 1;
        
        for quadrantIndex = 1:4
            %get quadrant image and indices
            [xOffset, yOffset, imageQuadrantOriginal] = ...
                CCDsegment4_GetQuadrant(dim1, dim2, imageOriginal, quadrantIndex);
            
            %smooth, if desired
            if opts.useSmoothing
                imageQuadrant = CCDsegment4_Smooth(opts, imageQuadrantOriginal);
            else
                imageQuadrant = imageQuadrantOriginal;
            end
            
            if UseBWLabel
                [T, E, ind, maskImage] = ...
                    CCDsegment4_bwlabel(imageQuadrant, xOffset, yOffset, threshold, ...
                    opts, T, E, ind);
            else
                [T, E, ind, maskImage] = ...
                    CCDsegment4_bwconncomp(imageQuadrant, xOffset, yOffset, threshold, ...
                    opts, T, E, ind);
            end
        end
        
        %cut off any extra cells or NaN's from T and E.
        T = T(1:ind-1);
        E = E(1:ind-1);
    end
    
else
    %segment entire image at once
    
    %initialize
    E = [];
    T = cell(0,1);
    ind = 1;
    
    %smooth, if desired
    if opts.useSmoothing
        imageFull = CCDsegment4_Smooth(opts, imageOriginal);
    else
        imageFull = imageOriginal;
    end
    
    if opts.useBWLabel
        [T, E, ind, maskImage] = ...
            CCDsegment4_bwlabel(imageFull, imageOriginal, 1, 1, threshold, ...
            opts, T, E, ind);
    else
        [T, E, ind, maskImage] = ...
            CCDsegment4_bwconncomp(imageFull, imageOriginal, 1, 1, threshold, ...
            opts, T, E, ind);
    end
    
    %cut off any extra cells or NaN's from T and E.
    T = T(1:ind-1);
    E = E(1:ind-1);
    
end

%{
    %~~~~~
    %create binary image. This is a logical array, where "true" marks
    %   pixels above the threshold, and "false" marks pixels below.
    lg = imageFull > threshold;
    
    %if we want to output the mask, we do it here...
    if maskFlag
        maskImage = bwmorph(lg,'dilate');    %mark neighbors too
        %this makes the assumption that neighborflag is true.
    end

    %label regions (tracks).
    %bwlabel marks contiguous groupings of "true" values with a unique number.
    %   L is an image with numbers indicating the segmented groupings;
    %   num is the total number of groups.
    [imageRegionsMask,numRegions] = bwlabel(lg,8);
    
    if numRegions > maxSegments
        %abort
        T = cell(0,1);
        E = [];
        return
    end
    
    T = cell(1,numRegions);
    E = nan(1,numRegions);
    
    for j=1:numRegions     %each track
        if sum(imageRegionsMask(:)==j) < NumPixelsThreshold     %there are not enough pixels in this grouping.
            continue    %skip to the next group
        end
        %otherwise, get a list of coordinates of the elements in this segmented group.
        [regionRows,regionCols] = find(imageRegionsMask==j);
        
        if edgeVeto && any(regionRows==1 | regionCols==1 | regionRows==dim1 | regionCols==dim2)
            %if any elements of this group are at the edge of the image,
            %   apply the edge veto if active.
            continue    %skip to the next group
        elseif any(regionRows==1 | regionCols==1 | regionRows==dim1 | regionCols==dim2)
            %if any elements of this group are at the edge of the image,
            %   make a note of it with edgeflag.
            edgeflag = true;
        else
            %not at the edge.
            edgeflag = false;
        end
        
        %mark neighbors if neighborflag is true.
        if NeighborLayers > 0
            regionMaskFull = bwmorph(imageRegionsMask==j,'dilate');   %save pixels above threshold, as well as neighboring pixels (8-conn)
        else
            regionMaskFull = imageRegionsMask==j;    %no neighbors
        end
        
        %get a list of coordinates of the elements in this group (now, including neighbors)
        [regionRows,regionCols] = find(regionMaskFull);
        
        %build the track image: it is padded with zeros outside the track.
        T{m}.img = zeros(max(regionRows)-min(regionRows)+1,max(regionCols)-min(regionCols)+1);  %create appropriately sized array
        
        regionMaskTrimmed = regionMaskFull(min(regionRows):max(regionRows),min(regionCols):max(regionCols));
        %"marksub" shows the pixels of the track in a rectangle the size of the track image
        T{m}.img(regionMaskTrimmed) = imageFull(regionMaskFull);     %copy marked pixels
        
        %record the other information of the track.
        T{m}.E = sum(T{m}.img(:));  %save energy of this track
        T{m}.x = min(regionRows)-1;            %rows are x dimension (vertical up)
        T{m}.y = min(regionCols)-1;            %columns are y dimension (horizontal away from collimator)
        %-1 compensates for being the offset (not the first row/column)
%         T{m}.q = NaN;             %which quadrant
        T{m}.edgeflag = edgeflag;   %whether this track is at the edge
                                    %   (if edgeveto is true, then this will only be false)
        E(m) = T{m}.E;  %in postcal mode, all energies are saved in the same list
        
        if ~isempty(subArrayParams)
            %subarray readout
            %adjust x and y to the coordinate system of the entire device, instead of the subarray.
            T{m}.y = T{m}.y + (subArrayParams(3)-subArrayParams(1)/2);  %columns
            T{m}.x = T{m}.x + (subArrayParams(4)-subArrayParams(2)/2);  %rows
        end
        m = m+1;
    end
    
    findNan = find(isnan(E),1);
    if ~isempty(findNan)
        %delete skipped ones
        Tnew = cell(1,findNan-1);
        for k=1:findNan-1
            Tnew{k} = T{k};
        end
        T = Tnew;
        E = E(1:findNan-1);
    end
end
%}
