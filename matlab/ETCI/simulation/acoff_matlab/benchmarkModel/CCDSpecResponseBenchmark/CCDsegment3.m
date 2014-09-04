function [T,E,mask] = CCDsegment3(varargin)
% function [T,E] = CCDsegment3(img,th)
% function [T,E] = CCDsegment3(img,th,'Mode',modestr)
% function [T,E] = CCDsegment3(img,th,'PixelThreshold',pixth)
% function [T,E] = CCDsegment3(img,th,'MirrorFlag',mirrorflag)
% function [T,E] = CCDsegment3(img,th,'NeighborFlag',NNflag)
% function [T,E] = CCDsegment3(img,th,'EdgeVeto',edgeveto)
% function [T,E] = CCDsegment3(img,th,'Subarray',[w,h,colc,rowc])
% function [T,E,mask] = CCDsegment3(...,'MaskFlag',maskflag)
%
%Inputs:
% img: input from CCDcal_median.
%       Black level should be normalized to 0.
% th: energy threshold for segmentation. Units defined by mode.
%     Precal mode: number of sigma of black level noise.
%     Postcal mode: keV
%      -March 10, 2010: trying 0.55 keV. (roughly equal to detector units 70)
%
% modestr: either 'precal' or 'postcal'
% pixth: number indicating minimum number of pixels in a segmented track.
% mirrorflag: logical; true means that mirrored flags will be removed.
% NeighborFlag: logical; true means that neighbor pixels will be included in track.
% edgeveto: logical; true means that tracks which touch the edge of the image will be vetoed.
% subarray: parameters to compensate for subarray readout of the CCD. 
%     w = box width, h = box height, colc = box column center, rowc = box rows center
%     bias width assumed to be zero
% 
% 
% Default options:
%   Mode = 'Postcal'
%   PixelThreshold = 4
%   MirrorFlag = false
%   NeighborFlag = true
%   EdgeVeto = false
%   MaskFlag = false
%
%Outputs:
% T: cell array of structures.
%   T{i}.img: image of track i.
%   T{i}.E: energy = sum(T{i}.img(:))
%   T{i}.x: the first row in the image, img(1,:), 
%               corresponds to img(x,:)
%   T{i}.y: the first column in the image, img(:,1), 
%               corresponds to img(:,y)
%           (x is vertical up in real space, 
%               y is horizontal away from collimator)
%   T{i}.edgeflag: marks if the track was neighboring the 
%       edge regions blacked out by CCDcal.
% E: list of energies. Either by quadrant (precal mode) or for the entire image.
% mask: binary mask of track pixels & neighbors.

%Input handling
if nargin < 2
    error('Need at least two input arguments.')
else
    img = varargin{1};
    th = varargin{2};
    if nargin==3
        warning('Three-argument function call deprecated... proceeding with img, th, pixth.')
        pixth = varargin{3};
    elseif mod(nargin,2)~=0
        error('Odd number of input arguments...')
    else
        %defaults
        precalflag = false;
        pixth = 4;
        mirrorflag = false;
        NNflag = true;
        quadrantstraddleveto = false;    %i.e. don't veto tracks crossing quadrant edges, by default
        edgeveto = true;
        subparflag = false;
        maskflag = false;
        
        for i = 3:2:nargin-1
            curarg = varargin{i+1};
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
            end
        end
    end
end

if maskflag && edgeveto
    error('Maskflag and Edgeveto should not be used together...')
end

T = [];

% d1 = 1454;
% d2 = 726;
[d1,d2] = size(img);
% th = 70;

m = 1;  %count saved tracks only

if quadrantstraddleveto
    if mirrorflag
        mmask = false(d1/2,d2/2);
    end
    if precalflag
        E = cell(1,4);
        thsig = th;
    else
        E = [];
    end
    
    %One quadrant at a time.
    for q=1:4
        switch q
            %clockwise order from upper left
            case 4
                xi = 1:d1/2;
                yi = 1:d2/2;
            case 1
                xi = d1/2+1:d1;
                yi = 1:d2/2;
            case 3
                xi = 1:d1/2;
                yi = d2/2+1:d2;
            case 2
                xi = d1/2+1:d1;
                yi = d2/2+1:d2;
        end
        imgq = img(xi,yi);
        
        if precalflag
            %fit background to determine threshold
            med1 = median(imgq(:));
            std1 = std(imgq(:));
            [n,x] = hist(imgq(:), (med1-6*std1):(std1/10):(med1+10*std1));
            opts1 = fitoptions('Method','LinearLeastSquares',...
                'Upper',[1.5*max(n),med1+1,sqrt(2)*std1*1.2],...
                'Lower',[0.8*max(n),med1-1,sqrt(2)*std1*0.8]);
            fit1 = fit(x',n','gauss1',opts1);
            th = med1 + thsig * fit1.c1 / sqrt(2);
        end
        %create binary image
        lg = imgq > th;
        if maskflag
            mask(xi,yi) = bwmorph(lg,'dilate');
        end
        
        %label regions (tracks)
        [L,num] = bwlabel(lg,8);

        m2 = 1; %energy index

        for j=1:num     %each track
            if sum(L(:)==j) < pixth ...     %below pixel threshold
                    || (mirrorflag && any(L(:)==j & mmask(:)))  %or vetoed by mirror
                continue
            end
            [r1,c1] = find(L==j);
            if any(r1==1 | c1==1 | r1==size(imgq,1) | c1==size(imgq,2)) %regardless of edgeveto setting
                continue
            end

            if NNflag
                mark = bwmorph(L==j,'dilate');   %save pixels above threshold, as well as neighboring pixels (8-conn)
    %             if any(mark(:) & imgq(:)==0)   %img==0 corresponds to edge regions which are blacked out by CCDcal
    %                 %a pixel in the edge region was marked (as a neighbor)
    %                 mark = mark & imgq>0;
    %                 edgeflag(m) = true;
    %             else
    %                 edgeflag(m) = false;
    %             end
                
            else
                mark = L==j;    %no neighbors
            end

            [r,c] = find(mark);

            T{m}.img = zeros(max(r)-min(r)+1,max(c)-min(c)+1);  %create appropriately sized array
            marksub = mark(min(r):max(r),min(c):max(c));
            T{m}.img(marksub) = imgq(mark);     %copy marked pixels
            T{m}.E = sum(T{m}.img(:));  %save energy of this track
            T{m}.x = min(r)+xi(1)-1;    %rows are x dimension (vertical up)
            T{m}.y = min(c)+yi(1)-1;    %columns are y dimension (horizontal away from collimator)
            T{m}.q = q;             %which quadrant
            T{m}.edgeflag = false;      %all edge tracks have already been vetoed
            if precalflag
                E{q}(m2) = T{m}.E;
            else
                E(m) = T{m}.E;
            end
            
            if subparflag
                %subarray readout
                %adjust x and y
                T{m}.y = T{m}.y + (subpar(3)-subpar(1)/2);  %columns
                T{m}.x = T{m}.x + (subpar(4)-subpar(2)/2);  %rows
            end
            m = m+1;
            m2 = m2+1;
        end
        
        if mirrorflag
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
        
    end
    
else
    %segment entire image at once
    
    E = [];
    %create binary image
    lg = img > th;
    
    if maskflag
        mask = bwmorph(lg,'dilate');
    end

    %label regions (tracks)
    [L,num] = bwlabel(lg,8);
    
    for j=1:num     %each track
        if sum(L(:)==j) < pixth     %below pixel threshold
            continue
        end
        [r1,c1] = find(L==j);
        if edgeveto && any(r1==1 | c1==1 | r1==d1 | c1==d2)
            continue
        elseif any(r1==1 | c1==1 | r1==d1 | c1==d2)
            edgeflag = true;
        else
            edgeflag = false;
        end
        
        if NNflag
            mark = bwmorph(L==j,'dilate');   %save pixels above threshold, as well as neighboring pixels (8-conn)
%             if any(mark(:) & imgq(:)==0)   %img==0 corresponds to edge regions which are blacked out by CCDcal
%                 %a pixel in the edge region was marked (as a neighbor)
%                 mark = mark & imgq>0;
%                 edgeflag(m) = true;
%             else
%                 edgeflag(m) = false;
%             end
        else
            mark = L==j;    %no neighbors
        end
        
        [r,c] = find(mark);
        
        T{m}.img = zeros(max(r)-min(r)+1,max(c)-min(c)+1);  %create appropriately sized array
        marksub = mark(min(r):max(r),min(c):max(c));
        T{m}.img(marksub) = img(mark);     %copy marked pixels
        T{m}.E = sum(T{m}.img(:));  %save energy of this track
        T{m}.x = min(r);            %rows are x dimension (vertical up)
        T{m}.y = min(c);            %columns are y dimension (horizontal away from collimator)
%         T{m}.q = NaN;             %which quadrant
        T{m}.edgeflag = edgeflag;
        E(m) = T{m}.E;
        
        if subparflag
            %subarray readout
            %adjust x and y
            T{m}.y = T{m}.y + (subpar(3)-subpar(1)/2);  %columns
            T{m}.x = T{m}.x + (subpar(4)-subpar(2)/2);  %rows
        end
        m = m+1;
    end
    
end
