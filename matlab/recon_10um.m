dataDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/10um';
cd(dataDir);
fileList = dir('./*.fits');
fileNum = length(fileList);
for i = 1:fileNum
    track = fitsread(fileList(i).name);
    recon_tmp = HybridTrack(track);
    if ~isfield(recon_tmp, 'err')
        recon = struct;
        recon.alpha = recon_tmp.alpha;
        recon.beta = recon_tmp.beta;
        % convert relative start position to absolute position
        h = fitsinfo(fileList(i).name);
        lowerLeft_rel = [h.PrimaryData.Keywords{14}, h.PrimaryData.Keywords{15}];
        startPix = recon_tmp.EdgeSegments.startCoordinatesPix;
        recon.startPix = lowerLeft_rel + startPix;
        % save to file, then delete
        saveName = sprintf('%d.mat', i);
        save(saveName, 'recon', '-v7.3');
        if mod(i, 500) == 0
            i
        end
        clear recon, recon_tmp, track, h;
    else
        clear recon_tmp, track;
    end
end

