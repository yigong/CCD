function [T,E,mask,noise] = CCDsegment3(varargin)
% function [T,E] = CCDsegment3(img,th)
% function [T,E] = CCDsegment3(img,th,'Mode',modestr)
% function [T,E] = CCDsegment3(img,th,'PixelThreshold',pixth)
% function [T,E] = CCDsegment3(img,th,'NeighborFlag',NNflag)
% function [T,E] = CCDsegment3(img,th,'EdgeVeto',edgeveto)
% function [T,E] = CCDsegment3(img,th,'Subarray',[w,h,colc,rowc])
% function [T,E,mask] = CCDsegment3(...,'MaskFlag',maskflag)
%
%Requires the Image Processing toolbox.
%
%Inputs:
% img: input image array; should come from CCDcal_median or similar.
%   The black level should be normalized to 0. In 'precal' mode, the image
%   should be in detector units (gain = 1); in 'postcal' mode, the image should
%   be in keV.
% th: energy threshold for segmentation. Units defined by mode.
%     Precal mode: units are [sigma of black level noise]; Yigong found that 4 is best.
%     Postcal mode: units are keV. Yigong found that 0.0715 is best for Am-241 data,
%       and higher-energy tracks might suggest an even lower threshold (haven't tried).
%
% modestr: either 'precal' or 'postcal'
%   'Precal' indicates that the gain of the device was unknown, so the images
%   are in detector units.
%   'Postcal' (default) means that the image is calibrated in keV.
% pixth: number indicating minimum number of pixels allowed in a segmented track.
%   If a track is found with n_pixels < pixth, it will not be stored.
% mirrorflag: logical; true means that mirrored flags will be removed.
%   **this is obsolete, since we fixed the wiring issue which caused mirrored tracks**
% NeighborFlag: logical; true means that neighbor pixels will be included in track.
%   Default is true. NeighborFlag = True has been found to improve energy resolution
%   in the Am-241 photopeak.
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
% 
% Default options:
%   Mode = 'Postcal'
%   PixelThreshold = 4
%   MirrorFlag = false
%   NeighborFlag = true
%   EdgeVeto = true
%   MaskFlag = false
%   maxSegments = +Inf
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

%% Input handling. 
%After the first two required arguments (img, th), look for parameter names 
%   and then read the parameter value from the next argument.
%There are some smarter ways to do this, but as it is, each parameter name
%   must be paired with the parameter value.
if nargin < 2
    error('Need at least two input arguments.')
else
    img = varargin{1};
    th = varargin{2};
    if mod(nargin,2)~=0
        %every optional parameter name should be paired with the parameter value.
        %if the number of arguments is odd, then this is not the case.
        error('Odd number of input arguments...')
    else
        %Set default values of optional parameters.
        %If they are specified in the argument list, then we will overwrite them.
        precalflag = false; %precalflag is true iff mode is 'precal'
        pixth = 4;
        mirrorflag = false; %mirrorflag is obsolete now.
        NNflag = true;      %neighborflag
        quadrantstraddleveto = false;    %i.e. don't veto tracks crossing quadrant edges, by default
        %I need to take a second look at QuadrantStraddleVeto . . I'm not sure it does anything right now.
        % QuadrantStraddleVeto is set as a result of other parameters in the inputs.
        %  #futurework
        edgeveto = true;    %also need to look at this again, i think it is written for an old CCDcal code. #futurework
        subparflag = false; %use subarray parameters
        maskflag = false;   %output mask array
        maxSegments = +Inf;
        
        for i = 3:2:nargin-1    %start from 3 because first two parameters are required
                    %read every second input argument, as a parameter name
            curarg = varargin{i+1}; %this is the value associated with the current parameter
            %look at the parameter name, identify it, and make sure that the 
            %   argument value is the correct data type.
            switch lower(varargin{i})
                case 'mode'
                    if strcmpi(curarg,'precal')
                        precalflag = true;
                        quadrantstraddleveto = true;    %veto tracks crossing quadrant edges, if we are performing a calibration
                    elseif strcmpi(curarg,'postcal')
                        precalflag = false;
                    else
                        error('Could not understand given mode. Mode should be ''precal'' or ''postcal'' ')
                    end
                case 'pixelthreshold'
                    if isnumeric(curarg) && length(curarg)==1
                        pixth = curarg;
                    else
                        error('Could not understand the given pixel threshold.')
                    end
                case 'mirrorflag'
                    if islogical(curarg) && length(curarg)==1
                        mirrorflag = curarg;
                    elseif isnumeric(curarg) && length(curarg)==1
                        mirrorflag = logical(curarg);
                    else
                        error('Could not understand the given mirror flag.')
                    end
                    if mirrorflag
                        quadrantstraddleveto = true;    %veto tracks crossing quadrant edges, if we are performing a calibration
                    end
                case 'neighborflag'
                    if islogical(curarg) && length(curarg)==1
                        NNflag = curarg;
                    elseif isnumeric(curarg) && length(curarg)==1
                        NNflag = logical(curarg);
                    else
                        error('Could not understand the given neighbor flag.')
                    end
                case 'edgeveto'
                    if islogical(curarg) && length(curarg)==1
                        edgeveto = curarg;
                    elseif isnumeric(curarg) && length(curarg)==1
                        edgeveto = logical(curarg);
                    else
                        error('Could not understand the given edge veto.')
                    end
                case 'subarray'
                    s = size(curarg);
                    if isnumeric(curarg) && all(s==[1,4])
                        subpar = curarg;
                        subparflag = true;
                    elseif isnumeric(curarg) && all(s==[4,1])   %column format?
                        subpar = curarg';
                        subparflag = true;
                    elseif isempty(curarg)  %no subarray
                        subparflag = false;
                    else
                        error('Could not understand the subarray parameters.')
                    end
                case 'maskflag'
                    maskflag = true;
                case 'maxsegments'
                    maxSegments = curarg;
                    
            end
        end
    end
end

if maskflag && edgeveto
    error('Maskflag and Edgeveto should not be used together...')
    %why? I don't know. These flags are old and outdated I think. #futurework
end

%% 
mask = false(size(img));
%measure dimensions from the CCD image given
[d1,d2] = size(img);

m = 1;  %this is an index for tracks that are saved to the "T" cell array

if quadrantstraddleveto
    %Segment each quadrant separately.
    if mirrorflag
        mmask = false(d1/2,d2/2);
    end
    if precalflag
        %list energy values separately by quadrant, not all in one list.
        %because each quadrant will have a slightly different gain value.
        E = cell(1,4);
        thsig = th;
        noise = zeros(1,4);
    else
        %energy values from all quadrants will all be in one list
        E = [];
    end
    T_temp = cell(1,4);
    for q=1:4   %each quadrant
        T = cell(1,100);

        switch q
            %Quadrants are numbered in clockwise order from upper left
            case 4
                %bottom left
                xi = 1:d1/2;
                yi = 1:d2/2;
            case 1
                %top left
                xi = d1/2+1:d1;
                yi = 1:d2/2;
            case 3
                %bottom right
                xi = 1:d1/2;
                yi = d2/2+1:d2;
            case 2
                %top right
                xi = d1/2+1:d1;
                yi = d2/2+1:d2;
        end
        %grab the image of one quadrant
        imgq = img(xi,yi);

            %the threshold is given relative to the noise (in sigma) of the black level.
            %we histogram all the pixels, which should be completely dominated by
            %  black level instead of tracks, and fit this histogram to get the 
            %  black level noise width. That tells us the numerical value of the 
            %  threshold.
        if precalflag
            %%%%%%%%%%%%%%%%%%% Updated by Yigong 07/18/2013
            noiseLowerLimit = -40; % in ADC units
            noiseUpperLimit = 40; % in ADC units
            [noiseCentroid, noiseStd] = noiseFitSeg(imgq(:), noiseLowerLimit, noiseUpperLimit);
            noise(q) = noiseStd;
            th = noiseCentroid + thsig * noiseStd;
            %%%%%%%%%%%%%%%%%%%
        end
        %create binary image. This is a logical array, where "true" marks
        %   pixels above the threshold, and "false" marks pixels below.
        lg = imgq > th;
        %if we want to output the mask, we do it here...
        if maskflag
            mask(xi,yi) = bwmorph(lg,'dilate'); %mark neighbors too
            %this makes the assumption that neighborflag is true.
        end
        %label regions (tracks).
        %bwlabel marks contiguous groupings of "true" values with a unique number.
        %   L is an image with numbers indicating the segmented groupings;
        %   num is the total number of groups.
        [L,num] = bwlabel(lg,8);
        m = 1;
        m2 = 1; %index for energy variable
        
        if num > maxSegments
            %abort
            T = cell(0,1);
            E = [];
            return
        end
        T = cell(1,num);
        E{q} = nan(1,num);
    
        for j=1:num     %each track
            if sum(L(:)==j) < pixth ...     %there are not enough pixels in this grouping.
                    || (mirrorflag && any(L(:)==j & mmask(:)))  %or vetoed by mirror (obsolete)
                continue    %skip to the next group
            end
            %otherwise, get a list of coordinates of the elements in this segmented group.
            [r1,c1] = find(L==j);
            if any(r1==1 | c1==1 | r1==size(imgq,1) | c1==size(imgq,2)) %regardless of edgeveto setting, in this mode.
                continue    %if any elements of this group are at the edge of the image,
                            %   skip to the next group. 
            end

            if NNflag
                %we need to mark the neighboring pixels as part of the same track.
                mark = bwmorph(L==j,'dilate');   %save pixels above threshold, as well as neighboring pixels (8-conn)
            else
                mark = L==j;    %no neighbors
            end
            %"mark" indicates the pixels out of the entire CCD image which belong to this grouping.
            
            %get a list of coordinates of the elements in this group (now, including neighbors)
            [r,c] = find(mark);
            
            %build the track image: it is padded with zeros outside the track.
            T{m}.img = zeros(max(r)-min(r)+1,max(c)-min(c)+1);  %create appropriately sized array
            
            marksub = mark(min(r):max(r),min(c):max(c));    
            %"marksub" shows the pixels of the track in a rectangle the size of the track image
            T{m}.img(marksub) = imgq(mark);     %copy the marked pixels
            
            %record the other information of the track.
            T{m}.E = sum(T{m}.img(:));  %save energy of this track
            T{m}.x = min(r)+xi(1)-2;    %rows are x dimension (vertical up)
            T{m}.y = min(c)+yi(1)-2;    %columns are y dimension (horizontal away from collimator)
            %-2 compensates both for the indexing of xi,yi as well as being the offset (not the first row/column)
            T{m}.q = q;             %which quadrant
            T{m}.edgeflag = false;      %all edge tracks have already been vetoed in this mode.
            if precalflag
                E{q}(m2) = T{m}.E;  %in precal mode, save energies separately from each quadrant
            else
                E(m) = T{m}.E;      %in postcal mode, all energies are saved in the same list
            end
            
            if subparflag
                %subarray readout is specified.
                %adjust x and y to the coordinate system of the entire device, instead of the subarray.
                T{m}.y = T{m}.y + (subpar(3)-subpar(1)/2);  %columns
                T{m}.x = T{m}.x + (subpar(4)-subpar(2)/2);  %rows
            end
            m = m+1;
            m2 = m2+1;
        end
        
        if mirrorflag
            %mirror the mask pattern over to the next quadrant (obsolete)
            switch q
                case 1
                    mmask = lg(:,end:-1:1);
                case 2
                    mmask = lg(end:-1:1,:);
                case 3
                    mmask = lg(:,end:-1:1);
                case 4
                    mmask = false(size(lg));
            end
        end
        
        findNan = find(isnan(E{q}),1);
        if ~isempty(findNan)
            %delete skipped ones
            Tnew = cell(1,findNan-1);
            for k=1:findNan-1
                Tnew{k} = T{k};
            end
            T_temp{q} = Tnew;
            E{q} = E{q}(1:findNan-1);
        else
            T_temp{q} = {};
        end
    end
    clear T
    T = [T_temp{1},T_temp{2},T_temp{3},T_temp{4}];

else
    %segment entire image at once
    
    %create binary image. This is a logical array, where "true" marks
    %   pixels above the threshold, and "false" marks pixels below.
    lg = img > th;
    
    %if we want to output the mask, we do it here...
    if maskflag
        mask = bwmorph(lg,'dilate');    %mark neighbors too
        %this makes the assumption that neighborflag is true.
    end

    %label regions (tracks).
    %bwlabel marks contiguous groupings of "true" values with a unique number.
    %   L is an image with numbers indicating the segmented groupings;
    %   num is the total number of groups.
    [L,num] = bwlabel(lg,8);
    
    if num > maxSegments
        %abort
        T = cell(0,1);
        E = [];
        return
    end
    
    T = cell(1,num);
    E = nan(1,num);
    
    for j=1:num     %each track
        if sum(L(:)==j) < pixth     %there are not enough pixels in this grouping.
            continue    %skip to the next group
        end
        %otherwise, get a list of coordinates of the elements in this segmented group.
        [r1,c1] = find(L==j);
        
        if edgeveto && any(r1==1 | c1==1 | r1==d1 | c1==d2)
            %if any elements of this group are at the edge of the image,
            %   apply the edge veto if active.
            continue    %skip to the next group
        elseif any(r1==1 | c1==1 | r1==d1 | c1==d2)
            %if any elements of this group are at the edge of the image,
            %   make a note of it with edgeflag.
            edgeflag = true;
        else
            %not at the edge.
            edgeflag = false;
        end
        
        %mark neighbors if neighborflag is true.
        if NNflag
            mark = bwmorph(L==j,'dilate');   %save pixels above threshold, as well as neighboring pixels (8-conn)
        else
            mark = L==j;    %no neighbors
        end
        
        %get a list of coordinates of the elements in this group (now, including neighbors)
        [r,c] = find(mark);
        
        %build the track image: it is padded with zeros outside the track.
        T{m}.img = zeros(max(r)-min(r)+1,max(c)-min(c)+1);  %create appropriately sized array
        
        marksub = mark(min(r):max(r),min(c):max(c));
        %"marksub" shows the pixels of the track in a rectangle the size of the track image
        T{m}.img(marksub) = img(mark);     %copy marked pixels
        
        %record the other information of the track.
        T{m}.E = sum(T{m}.img(:));  %save energy of this track
        T{m}.x = min(r)-1;            %rows are x dimension (vertical up)
        T{m}.y = min(c)-1;            %columns are y dimension (horizontal away from collimator)
        %-1 compensates for being the offset (not the first row/column)
%         T{m}.q = NaN;             %which quadrant
        T{m}.edgeflag = edgeflag;   %whether this track is at the edge
                                    %   (if edgeveto is true, then this will only be false)
        E(m) = T{m}.E;  %in postcal mode, all energies are saved in the same list
        
        if subparflag
            %subarray readout
            %adjust x and y to the coordinate system of the entire device, instead of the subarray.
            T{m}.y = T{m}.y + (subpar(3)-subpar(1)/2);  %columns
            T{m}.x = T{m}.x + (subpar(4)-subpar(2)/2);  %rows
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
