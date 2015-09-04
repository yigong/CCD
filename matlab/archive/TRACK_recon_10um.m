matlabpool open 12
dataDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/10um';
cd(dataDir);
fileList = dir('./*.fits');
fileNum = length(fileList);
parfor i = 1:fileNum
    track = fitsread(fileList(i).name);
    recon_tmp = HybridTrack(track);
    ['recon # ', num2str(i)]
    if ~isfield(recon_tmp, 'err')
        ['save # ', num2str(i)]
        recon = zeros(1,4);
        recon(1) = recon_tmp.alpha;
        recon(2) = recon_tmp.beta;
        % convert relative start position to absolute position
        h = fitsinfo(fileList(i).name);
        minRow_track = h.PrimaryData.Keywords{14}
        minCol_track = h.PrimaryData.Keywords{15}
        startPix = recon_tmp.EdgeSegments.startCoordinatesPix;
        recon(3) = minRow_track + startPix(1)
        recon(4) = minCol_track + startPix(2)
       % % save to file, then delete
        saveName = sprintf('%d.mat', i);
        saveParFor(saveName, recon, '-v7.3');
        % free up memory
        recon = [];
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
