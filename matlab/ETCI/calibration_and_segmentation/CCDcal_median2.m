function CCDcal_median2(basedirin,dirlist,flist,filesuffix,basedirout,gain,replaceflag)
%function CCDcal_median2(basedirin,dirlist,flist,filesuffix,basedirout,gain,replaceflag)
%
% (v2: now each quadrant of each image is adjusted to its mean, after medianing)
%
%   Perform energy calibration on a series of CCD images.
%   The code is set up for images that may be split across a series of folders,
% e.g. the images from the coincidence experiment which are divided by hour.
% The images are arranged in chronological order, then batches of 100 images 
% at a time are calibrated using a median operation on each individual pixel
% to compensate for any black level variations across the image.
%   The median values are subtracted from each individual image, and then the 
% result multiplied by the gain values of each quadrant, to get a calibrated 
% image.
%
%   For black level normalization with unknown gain values, use gains of 1.
%   To disable the output of a quadrant, give a gain of 0.
%   Quadrants are numbered from top left (d1/2:d1,1:d2/2), clockwise, to 
% bottom left (1:d1/2,1:d2/2).
%
%
% inputs:
%
% basedirin: A string indicating the part of the file path that is
%   common to all input files. e.g. '/mnt/grimdata5/CCD/co60_20120730/'
% dirlist: A cell array of strings, each indicating a subpath within 
%   basedirin. **There is one cell for each .fits file, which contains the
%   subpath where that .fits file is located!**
% flist: File structure array, the result of 'dir' command in each subpath,
%   concatenated together. The path should be contained in 'basedirin' and
%   'dirlist'. 
%   Then the full path/name of input file 'i' is given by
%   fullfile(basedirin, dirlist{i}, flist(i).name)
%   (Each element of flist must include fields 'name' and 'datenum.')
% filesuffix: A string which will be appended to each output filename,
%   before the file extension. 
%   Thus, if filesuffix = '_suf', then 'img01.fits' becomes 'img01_suf.mat'.
% basedirout: A string indicating the location to put output files in.
%   The output files will be placed within basedirout, in subfolders described
%   by dirlist.
%   So the full path/name of output file 'i' is given by something like
%   fullfile(basedirout, dirlist{i}, [flist(i).name(1:n), filesuffix, '.mat'])
%   where n is the filename position before '.fits'
% gain: A 1x4 vector of gain values, for quadrants of CCD from top left clockwise.
%   Use 1 if the gain is not calibrated yet, in order to get tracks in
%   detector units. Use 0 to disable a quadrant.
% replaceflag: flag for replacing images already calibrated. Otherwise, will skip these.
% 
% added comments, Aug 31, 2012

%use a set of this many images to get the black level baseline for each pixel.
%(This is a minimum)
batchsize = 100;

% median_adjust = true;  %true: adjust each individual image so that the median *of each quadrant* is at 0.
% ...this is always true now.

%load the first file to measure the dimensions of the image
f = fitsread(fullfile(basedirin,dirlist{1},flist(1).name));
[d1,d2] = size(f);
%d1 is height (in pixels)
%d2 is width (in pixels)
disp(['Found images of size ',num2str(d1),' by ',num2str(d2)])

%start timer. t1 is the start value in seconds.
try
    t1 = toc;   
    %if toc returns a value, that means that tic is already running, 
    %so just measure the current stopwatch time instead of starting a new one
catch
    %tic is not running, so we should start it.
    tic;
    t1 = 0;
end

%check which files we really need to calibrate.
%at the same time, make note of the modification timestamp on each file,
% so we can make sure we calibrate them in chronological order.
% (normally, the files are numerically named in order, so flist returns
%  a list that is in chronological order already.)
flg = false(size(flist));   
%flg is a logical value for each file in flist.
%   if flg(i) is true, it indicates the file can be skipped. however, the
%   skipped files might still be used for the median operation to measure
%   the black level for nearby non-skipped files. see "while" loop below.
%   (default is not to skip anything.)
for i=1:length(flist)
    ts(i) = flist(i).datenum;   %note the timestamp
    if ~replaceflag && ...      %if replaceflag is false, then all files will be calibrated regardless.
            ~isempty(dir(fullfile(basedirout,dirlist{i},convertfilename(flist(i).name,filesuffix))))
        %output mat file already exists
        flg(i) = true;  %so we will skip this file
    end
end

%sort file list by timestamp
[ts,ind] = sort(ts);
flist = flist(ind); %now flist is in timestamp order.

%rearrange dirlist to match the sorted flist.
%(yeah, this dirlist thing should be optional rather than hardcoded.... #futurework)
tmp = dirlist;  %we can't reference dirlist{ind} because it is a cell array.
%so we need a buffer cell array as we rearrange things; 
% tmp is the buffer that we read dirlist entries from.
for i = 1:length(dirlist)
    dirlist{i} = tmp{ind(i)};
end

%rearrange the flg (skip flags) as well.
flg = flg(ind);

%now, set up batches of batchsize for median-ing. a batch should be batchsize
%   number of files unless there is an uneven number at the end, in which case
%   the batch will be more than batchsize.
%   In other words, batchsize gives the minimum size of the batch.

%batch index variable "b" will list the index of the first file of the batch in column 1,
%   and the index of the last file of the batch in column 2.
if length(ts)/batchsize >= 2    %this means we have multiple batches.
    maxind = floor(length(ts)/batchsize);   %this is the number of batches there will be
    b = zeros(maxind,2);    %set up batch index variable "b"
    b(:,1) = 1:batchsize:batchsize*(maxind-1)+1;    %start index is always (N-1)*batchsize + 1
    b(:,2) = batchsize:batchsize:length(ts);        %finish index is usually N*batchsize
    b(size(b,1),2) = length(ts);                    %but last finish index is the end of the list.
else    %this means we only have one batch. incidentally, it is allowed to 
        % be smaller than batchsize, which could be a bad thing if
        % this code is misused.
    b(1,1) = 1; %start at index 1
    b(1,2) = length(ts);    %finish at last index
end

%if any batch can be skipped entirely, then remove the entry from the 
%   batch index variable "b".
i=1;
while i <= size(b,1)
    if all(flg(b(i,1):b(i,2)))
        %all files can be skipped
        b(i,:) = [];
    else
        i = i+1;
    end
end

%output info. (batch setup takes <= 1 second)
t2 = toc;
disp(['batch setup of ',num2str(length(ts)),' files took ',num2str(t2-t1),' seconds'])
t1 = t2;

%now, we calibrate each batch.
for i=1:size(b,1)   %each batch
    bvec = b(i,1):b(i,2);   %list of file indices in this batch
    
    %This 3D array stores every image in the batch. In other words, it is a
    % "stack" of ~batchsize number of CCD images.
    % I use uint16 instead of double. Uint16 is the format of the *.fits data,
    % and takes less space than a double. (This 3D array is ~100 MB of memory.)
    
    m = uint16(zeros(d1,d2,length(bvec))); %initialize the 3D array
    
    %load each image
    for j = 1:length(bvec)  %each image in batch
        %read fits file
        f = fitsread(fullfile(basedirin,dirlist{bvec(j)},flist(bvec(j)).name));
        %store as a uint16
        f = uint16(f);
        %check image dimensions - if they don't match then we have a problem.
        if size(f,1) ~= d1 || size(f,2) ~= d2
            error(['Bad image size: ',fullfile(basedirin,dirlist{bvec(j)},flist(bvec(j)).name)])
        end
        %otherwise, we copy this image into its place in the "stack".
        m(:,:,j) = f;
    end
    
    %output info. (loading files usually takes ~100 seconds)
    t2 = toc;
    disp(['Images ',num2str(b(i,1)),' to ',num2str(b(i,2)),' loaded in ',num2str(t2-t1),' seconds'])
    t1 = t2;
    
    %median the batch. The median operates in dimension 3, or the "height" of the "stack".
    %So each pixel is compared to the same position pixel in the other images,
    % and the result is a median image containing the median value of each pixel.
    med = median(m,3);
    
    %output info. (medianing the data usually takes ~30 seconds)
    t2 = toc;
    disp(['Images ',num2str(b(i,1)),' to ',num2str(b(i,2)),' medianed in ',num2str(t2-t1),' seconds'])
    t1 = t2;
    
    %finish calibrating each image, and save to disk.
    for j=1:length(bvec)
        %use the median image to adjust this individual image, and apply gain.
        img = calibrateCCDimage(m(:,:,j),med,d1,d2,gain); 
        %create the directory structure if needed.
        if isempty(dir(fullfile(basedirout,dirlist{bvec(j)})))
            mkdir(fullfile(basedirout,dirlist{bvec(j)}));
        end
        %write to disk as a *.mat file.
        save(fullfile(basedirout,dirlist{bvec(j)},convertfilename(flist(bvec(j)).name,filesuffix)),'img')
    end
    
    %output info. (saving files usually takes ~150 seconds)
    t2 = toc;
    disp(['Images ',num2str(b(i,1)),' to ',num2str(b(i,2)),' saved in ',num2str(t2-t1),' seconds'])
    t1 = t2;
    
    %get rid of some variables.
    clear m f img imgd med
end



function newname = convertfilename(filename,filesuffix)
    %change file extension from *.fits to *.mat
    newname = [filename(1:strfind(filename,'.')-1),filesuffix,'.mat'];

function img = calibrateCCDimage(img,med,d1,d2,gain)
    %this is where a single image is calibrated with pixel offsets and quadrant gains.
    
    %subtract black level from each pixel
    imgd = double(img) - double(med);
    
    img = imgd - median(imgd(:));   %this is no longer necessary since we are adjusting each quadrant.
    
    %quadrants are defined from top left, clockwise, as usual.
    q1 = imgd(d1/2+1:d1, 1:d2/2);
    q2 = imgd(d1/2+1:d1, d2/2+1:d2);
    q3 = imgd(1:d1/2, d2/2+1:1:d2);
    q4 = imgd(1:d1/2, 1:d2/2);
    %adjust median of each quadrant within the single image
    img(d1/2+1:d1, 1:d2/2) =    q1 - median(q1(:));
    img(d1/2+1:d1, d2/2+1:d2) = q2 - median(q2(:));
    img(1:d1/2, d2/2+1:1:d2) =  q3 - median(q3(:));
    img(1:d1/2, 1:d2/2) =       q4 - median(q4(:));
    %apply gain to each quadrant
    img(d1/2+1:d1, 1:d2/2) = img(d1/2+1:d1, 1:d2/2)         .* gain(1); %q1 = top left
    img(d1/2+1:d1, d2/2+1:d2) = img(d1/2+1:d1, d2/2+1:d2)   .* gain(2); %q2 = top right
    img(1:d1/2, d2/2+1:1:d2) = img(1:d1/2, d2/2+1:1:d2)     .* gain(3); %q3 = bottom right
    img(1:d1/2, 1:d2/2) = img(1:d1/2, 1:d2/2)               .* gain(4); %q4 = bottom left