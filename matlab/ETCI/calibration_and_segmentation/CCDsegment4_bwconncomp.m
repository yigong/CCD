function [T, E, ind, maskImage] = CCDsegment4_bwconncomp(imagePortion, imagePortionOriginal, xOffset, yOffset, threshold, opts, T, E, ind)
%function [T, E, ind, maskImage] = CCDsegment4_bwconncomp(imagePortion, imagePortionOriginal, xOffset, yOffset, threshold, opts, T, E, ind)
%
% Operated by CCDsegment4.

if length(threshold) == 2
    [dim1, dim2] = size(imagePortion);
    %create binary image.
    binImage_bottom = imagePortion(1:dim1/2, :) > threshold(1);
    binImage_top = imagePortion(dim1/2+1:end, :) > threshold(2);
    binImage = [binImage_bottom; binImage_top];
else 
    binImage = imagePortion > threshold;
end

%add layers of neighboring pixels, except the last layer.
%   this matches CCDsegment3 behavior, and prevents unnecessary restrictions
%   on larger pixel sizes, while keeping small-pixel, noisy tracks contiguous.
binImage = bwmorph(binImage,'dilate',max(1,opts.neighborLayers)); %in case neighborLayers==0.

%set up mask. tracks will be added one by one, after they are validated.
maskImage = false(size(binImage));

% %if we want to output the mask, we do it here...
% if maskFlag && NeighborsFlag
%     %add the last layer of neighbor pixels.
%     maskImage(xQuadrant,yQuadrant) = bwmorph(binImage,'dilate');
% end

%label regions (tracks).
connStruct = bwconncomp(binImage);

%check for maxSegments exceeded.
if connStruct.NumObjects > opts.maxSegments
    return
end

%filter out by pixel threshold.
abovePixelsThreshold = true(1,connStruct.NumObjects);
if opts.numPixelsThreshold > 0
    for i=1:connStruct.NumObjects
        abovePixelsThreshold(i) = length(connStruct.PixelIdxList{i}) >= opts.numPixelsThreshold;
    end
end

findAbovePixelsThreshold = find(abovePixelsThreshold);
if isempty(findAbovePixelsThreshold)
    %nothing to do.
    return
end

%initialize region mask variable, so we can just index it later.
regionMaskFull = false(connStruct.ImageSize);
%initialize outputs. They still might need trimming due to edgeVeto.
T{ind-1+length(findAbovePixelsThreshold)} = [];   %allocate up to this cell index of the cell array.
E(ind:ind-1+length(findAbovePixelsThreshold)) = nan;

%each track above pixels threshold. ("region")
for j=1:length(findAbovePixelsThreshold)
    %get a list of subscripts.
    [regionRows, regionCols] = ind2sub(connStruct.ImageSize, ...
        connStruct.PixelIdxList{findAbovePixelsThreshold(j)});
    
    %if any elements of this group are at the edge of the image, skip to the next group.
    if any(regionRows==1 | regionCols==1 | ...
            regionRows==connStruct.ImageSize(1) | regionCols==connStruct.ImageSize(2))
        if opts.edgeVeto
            continue
        else
            T{ind}.edgeflag = true;
        end
    end
    
    %clear the region mask
    regionMaskFull(:) = false;
    %mark the pixels of this region in a binary image
    regionMaskFull(connStruct.PixelIdxList{findAbovePixelsThreshold(j)}) = true;
    
    if opts.neighborLayers > 0
        %we need to mark the neighboring pixels as part of the same track.
        regionMaskFull = bwmorph(regionMaskFull,'dilate');   %save pixels above threshold, as well as neighboring pixels (8-conn)
    end
    %"regionMaskFull" indicates the pixels out of the entire CCD image which belong to this grouping.
    
    %now, add this track to the image mask.
    maskImage(regionMaskFull) = true;
%     maskImage = maskImage | regionMaskFull;   %same thing
    
    %get a list of coordinates of the elements in this group (now, including neighbors)
    [regionRows,regionCols] = find(regionMaskFull);
    
    %build the track image: it is padded with zeros outside the track.
    T{ind}.img = zeros(max(regionRows)-min(regionRows)+1,max(regionCols)-min(regionCols)+1);  %create appropriately sized array

    regionMaskTrimmed = regionMaskFull(min(regionRows):max(regionRows),min(regionCols):max(regionCols));    
    %"regionMaskTrimmed" shows the pixels of the track in a rectangle the size of the track image
    %   here is where we use original pixel data instead of the (possibly) smoothed data
    T{ind}.img(regionMaskTrimmed) = imagePortionOriginal(regionMaskFull);     %copy the marked pixels

    %record the other information of the track.
    T{ind}.E = sum(T{ind}.img(:));  %save energy of this track
    T{ind}.x = min(regionRows) + xOffset - 2;    %rows are x dimension (vertical up)
    T{ind}.y = min(regionCols) + yOffset - 2;    %columns are y dimension (horizontal away from collimator)
    %-2 compensates both for the indexing of xi,yi as well as being the offset (not the first row/column)
    
    %CCDsegment4: T{i}.q is not needed, because T is split into 4 cells, just like E.
%     T{ind}.q = quadrantIndex;             %which quadrant
    %edgeflag taken care of above, when edgeVeto is checked.
%     T{ind}.edgeflag = false;      %all edge tracks have already been vetoed in this mode.
    
    E(ind) = T{ind}.E;
    
    if ~isempty(opts.subArrayParams)
        %subarray readout is specified.
        %adjust x and y to the coordinate system of the entire device, instead of the subarray.
        T{ind}.y = T{ind}.y + (opts.subArrayParams(3) - opts.subArrayParams(1)/2);  %columns
        T{ind}.x = T{ind}.x + (opts.subArrayParams(4) - opts.subArrayParams(2)/2);  %rows
    end
    ind = ind+1;
end
