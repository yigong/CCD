function Ridge = RidgeFollow(img, Init, fPlot)
% compute alpha, beta using ridge following algorithm

% initialize 
Ntransect = 46;
width = zeros(Ntransect, 2);
isFinished = false;

rowNum = size(img, 1);
colNum = size(img, 2);
[colGV, rowGV] = meshgrid(1:colNum, 1:rowNum);
colGV = colGV + 0.5;
rowGV = rowGV + 0.5;
    
[rowRidge, colRidge] = Init.RFpos; % should be half-int start from 1.5
alphaG = Init.alphaGuess;
stepAngle = int16(Init.stepAngle);
col_0deg = -5:0.25:5;
row_0deg = zeros(size(col_0deg));
[lenGrid, cosGrid] = meshgrid(col_0deg, cosd([0:180])); % generate grid for rotation matrix
[lenGrid, sinGrid] = meshgrid(col_0deg, sind([0:180])); % generate grid for rotation matrix
rowsRotation = lenGrid .* sinGrid;
colsRotation = lenGrid .* cosGrid;

ridgePos = zeros(50, 2);
ridgePos(1, 1) = rowRidge;
ridgePos(1, 2) = colRidge;
for j = 2:50
    transectDeg = (stepAngle-90-Ntransect/2) : (stepAngle-90+Ntransect/2);
    transectDeg(transectDeg>180) = transectDeg(transectDeg>180) - 180;
    transectDeg(transectDeg<0)   = transectDeg(transectDeg<0) + 180;
    rowsTransect = rowsRotation(transectDeg, :) + rowRidge;
    colsTransect = colsRotation(transectDeg, :) + colRidge;
    transectEnrg = interp2(colGV, rowGV, img/16, colsTransect, rowsTransect);
    transectEnrg(isnan(transectEnrg)) = 0;
    for i = 1:size(transectEnrg, 1) 
        width(i,:) = fwhm(col_0deg, transectEnrg(i,:));
    end
    [~, iMinWidth] = min(width(:, 2));  % find minWidthCut index and the centroid index
    iCentroid = width(iMinWidth, 1);
    minWidthAngle = transectDeg(iMinWidth);
    alphaOpts = [minWidthAngle-90, minWidthAngle+90];
    [~, iAlpha] = min(abs(alphaOpts - alphaG));
    stepAngle = alphaOpts(iAlpha);
    rowRidge = rowsTransect(iMinWidth, iCentroid) + 0.25 * sind(double(stepAngle)); % rowsTransect(iMinWidth, iCentroid) gives centroid row 
    colRidge = colsTransect(iMinWidth, iCentroid) + 0.25 * cosd(double(stepAngle));
    ridgePos(j, :) = [rowRidge, colRidge];
    if transectEnrg(iMinWidth, iCentroid) < 0.5/16
        break;
    end
end 
pos = ridgePos(1:j, :);
Ridge.pos = pos;

if fPlot
load cmaps.mat;
cmap = cmaphotlog;
    PlotImage(img, cmap);
    hold on
    plot3(pos(:,2), pos(:,1), ones(size(pos(:,2))), 'markeredgecolor', 'g', 'markersize', 4)
end

end
