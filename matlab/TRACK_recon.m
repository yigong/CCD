function TRACK_recon(fitsDir, matFile)
    matlabpool open 12
    cd(fitsDir);
    fileList = dir('./*.fits');
    fileNum = length(fileList);
    result_arr = zeros(1, fileNum);
    parfor i = 1:fileNum
        fName = fileList(i).name;
        track = fitsread(fName);
        row_col = size(track);
        noise = randn(row_col) * 0.025;
        track = track + noise;
        h = fitsinfo(fName);
        dotIdx = strfind(fName, '.');
        trackID = str2num(fName(1:dotIdx-1));
        recon_tmp = HybridTrack(track, 'energyT', 1200.);
        ['recon # ', trackID]
        if ~isfield(recon_tmp, 'err')
            result_arr(i) = recon_tmp.alpha;
            recon_tmp = [];     % free up memory
            track = [];
            h = [];
       else
            recon_tmp = [];     % free up memory
            track = [];
            h = [];
        end
    end
    matlabpool close
    save(matFile, 'result_arr', '-v7')
return
