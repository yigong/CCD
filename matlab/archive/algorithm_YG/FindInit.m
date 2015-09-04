function Init = FindInit(img, pFlag)
% find the initial end

NeibrCt = ones(3); 
thr = 0.2;
% thr = Options.lowThresholdkeV;

binImg = +(img > thr);
thnImg = +bwmorph(binImg, 'thin', inf);
NeibrNO = conv2(thnImg, NeibrCt, 'same') - 1;
NeibrNO_thn = NeibrNO .* thnImg;
trackEnds = (NeibrNO_thn == 1); % binary image indicates ends locations
[rows, cols] = find(trackEnds); % find the row and col of ends
endsLinearIdc = find(trackEnds); % linear indices to extract Eends
psf = ones(6);
psf([1,6], [1,2,5,6]) = 0;
psf([2,5], [1,6]) = 0;
EnrgSumImg = conv2(img, psf, 'same');
endsEnrg = EnrgSumImg(endsLinearIdc)
[~, i] = min(endsEnrg);
r = rows(i);
c = cols(i);
endPos = [r+0.5, c+0.5];  % save the init row and col ||  r, c will change later
thnImg_copy = thnImg;
for iStep = 1:3
% for iStep = 1:Options.nStepBack
    thnImg_copy(r, c) = 0;
    neighbors = thnImg_copy(r-1:r+1, c-1:c+1);
    [r_rel, c_rel] = find(neighbors);
    r_rel = r_rel - 2;
    c_rel = c_rel - 2;
    r = r + r_rel;
    c = c + c_rel;
end
RFpos = [r+0.5, c+0.5];
stepAngle = (atan2(r_rel, c_rel) + pi)/pi*180; % in radius
alphaGuess = (atan2(RFpos(1)-endPos(1), RFpos(2)-endPos(2)) + pi)/pi*180; % in radius
% assign values to output
Init.endPos = endPos;
Init.RFpos = RFpos;
Init.stepAngle = stepAngle;
Init.alphaGuess = alphaGuess;

if pFlag
load cmaps.mat;
cmap = cmaphotlog;
    PlotImage(img, cmap);
    PlotImage(binImg, 'gray');
    PlotImage(thnImg, 'gray');
end







