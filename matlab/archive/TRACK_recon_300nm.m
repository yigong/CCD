matlabpool open 12
currentDir = pwd;
dataDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/300nm/fits_pixelP';
cd(dataDir);
fileList = dir('./*.fits');
fileNum = length(fileList);
parfor i = 1:fileNum
    fName = fileList(i).name;
    track = fitsread(fName);
    h = fitsinfo(fName);
    dotIdx = strfind(fName, '.');
    trackID = fName(1:dotIdx-1)
    eInit_true = h.PrimaryData.Keywords{11,2};
    recon_tmp = HybridTrack(track, 'energyT', eInit_true);
    ['recon # ', trackID]
    if ~isfield(recon_tmp, 'err')
        ['save # ', trackID]
        result = struct();
        % angle
        result.alphaM = recon_tmp.alpha;
        result.betaM = recon_tmp.beta;
        result.alphaT = h.PrimaryData.Keywords{8, 2};
        % position
        result.xT = h.PrimaryData.Keywords{9, 2};
        result.yT = h.PrimaryData.Keywords{10, 2};
        trackPosRow = h.PrimaryData.Keywords{6, 2};
        trackPosCol = h.PrimaryData.Keywords{7, 2};
        result.endRow = trackPosRow + recon_tmp.EdgeSegments.chosenEnd(1) - 1;
        result.endCol = trackPosCol + recon_tmp.EdgeSegments.chosenEnd(2) - 1;
        result.xM = result.endCol*10.5 - 37000/2; 
        result.yM = result.endRow*10.5 + 2000;
        % energy
        result.ET = eInit_true;
        result.EM = recon_tmp.energyM; 
        % save to file, then delete
        saveName = sprintf('./%s.mat', trackID);
        saveParFor(saveName, result, '-v7');
        % free up memory
        result = [];
        recon_tmp = [];
        track = [];
        h = [];
   else
        % free up memory
        recon_tmp = [];
        track = [];
        h = [];
    end
end
cd(currentDir);
matlabpool close
