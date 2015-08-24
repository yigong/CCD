function [T, E, ind, maskImage] = CCDsegment4_bwlabel(imagePortion, imagePortionOriginal, xOffset, yOffset, threshold, opts, T, E, ind)
%function [T, E, ind, maskImage] = CCDsegment4_bwlabel(imagePortion, imagePortionOriginal, xOffset, yOffset, threshold, opts, T, E, ind)
%
% Operated by CCDsegment4.

%create binary image.
binImage = imagePortion > threshold;
%add layers of neighboring pixels, except the last layer.
%   this matches CCDsegment3 behavior, and prevents unnecessary restrictions
%   on larger pixel sizes, while keeping small-pixel, noisy tracks contiguous.
binImage = bwmorph(binImage,'dilate',max(0,opts.neighborLayers-1)); %in case neighborLayers==0.

%set up mask. tracks will be added one by one, after they are validated.
maskImage = false(size(binImage));

% %if we want to output the mask, we do it here... 
% if maskFlag && NeighborsFlag
%     %add the last layer of neighbor pixels.
%     maskImage(xQuadrant,yQuadrant) = bwmorph(binImage,'dilate');
% end
%instead, add to the mask as each track passes the pixel threshold.

%label regions (tracks).

%bwlabel marks contiguous groupings of "true" values with a unique number.
%   imageRegionsMask is an image with numbers indicating the segmented groupings;
%   numRegions is the total number of groups.
[imageRegionsMask,numRegions] = bwlabel(binImage,8);

if numRegions > opts.maxSegments || numRegions==0
    %abort
    return
end
T{ind-1+numRegions} = [];   %allocate this portion of the array.
E(ind:ind-1+numRegions) = nan;

for j=1:numRegions     %each track
    if sum(imageRegionsMask(:)==j) < opts.numPixelsThreshold     %there are not enough pixels in this grouping.
        continue    %skip to the next group
    end
    %otherwise, get a list of coordinates of the elements in this segmented group.
    [regionRows,regionCols] = find(imageRegionsMask==j);
    if opts.edgeVeto && any(regionRows==1 | regionCols==1 | ...
            regionRows==size(imagePortion,1) | regionCols==size(imagePortion,2))
        continue    %if any elements of this group are at the edge of the image,
                    %   skip to the next group. 
    end
    
    if opts.neighborLayers > 0
        %we need to mark the neighboring pixels as part of the same track.
        regionMaskFull = bwmorph(imageRegionsMask==j,'dilate');   %save pixels above threshold, as well as neighboring pixels (8-conn)
    else
        regionMaskFull = imageRegionsMask==j;    %no neighbors
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
    T{ind}.edgeflag = false;      %all edge tracks have already been vetoed in this mode.
    
    E(ind) = T{ind}.E;
    
    if ~isempty(opts.subArrayParams)
        %subarray readout is specified.
        %adjust x and y to the coordinate system of the entire device, instead of the subarray.
        T{ind}.y = T{ind}.y + (opts.subArrayParams(3) - opts.subArrayParams(1)/2);  %columns
        T{ind}.x = T{ind}.x + (opts.subArrayParams(4) - opts.subArrayParams(2)/2);  %rows
    end
    ind = ind+1;
end
