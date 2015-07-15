matlabpool open 12
dataDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/1um/fits'
cd(dataDir);
fileList = dir('./*.fits');
fileNum = length(fileList);
parfor i = 1:fileNum
    track = fitsread(fileList(i).name);
    recon_tmp = HybridTrack(track);
    ['recon # ', num2str(i)]
    if ~isfield(recon_tmp, 'err')
        ['save # ', num2str(i)]
        var = zeros(1,8);
        var(1) = recon_tmp.alpha;
        var(2) = recon_tmp.beta;
        % convert relative start position to absolute position
        h = fitsinfo(fileList(i).name);
        minRow_track = h.PrimaryData.Keywords{6,2};
        minCol_track = h.PrimaryData.Keywords{7,2};
        startPix = recon_tmp.EdgeSegments.startCoordinatesPix;
        var(3) = minRow_track + startPix(1);
        var(4) = minCol_track + startPix(2);
        % save true values
        alphaTrue = h.PrimaryData.Keywords{8,2};
        xInit_true = h.PrimaryData.Keywords{9,2};
        yInit_true = h.PrimaryData.Keywords{10,2};
        eInit_true = h.PrimaryData.Keywords{11,2};
        var(5) = alphaTrue;
        var(6) = xInit_true;
        var(7) = yInit_true;
        var(8) = eInit_true;
       % % save to file, then delete
        saveName = sprintf('../mat/%d.mat', i);
        saveParFor(saveName, var, '-v7.3');
        % free up memory
        var= [];
        recon_tmp = [];
        track = [];
        h = [];
   else
        % free up memory
        recon_tmp = [];
        track = [];
    end
end
matlabpool close
