function algorithm_test(fitsDir, matDir)
    matlabpool open 12
    currentDir = pwd;
    cd(fitsDir);
    fileList = dir('*.fits');
    fileNum = length(fileList);
    alphaLR = zeros(size(fileList));
    alphaMD = zeros(size(fileList));
    EM = zeros(size(fileList));
    parfor i = 1:fileNum
        fName = fileList(i).name;
        dotIdx = strfind(fName, '.');
        trackID = int16(str2num(fName(1:dotIdx-1)));
        track = fitsread(fName);
        recon = HybridTrack(track, 'energyT', 500);
        if ~isfield(recon_tmp, 'err')
            alphaLR(trackID) = recon.alpha_LR;
            alphaMD(trackID) = recon.alpha;
            EM(trackID) = recon.energyM;
        end
    end
    saveFile = [matDir, '/algorithm_test.mat'];
    save(saveFile, 'alphaLR', 'alphaMD', 'EM');


