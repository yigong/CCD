dataDir = '/global/scratch/ygzhang/G4BeamDiag/step_size/out_LRC/3um/mat';
cd(dataDir);
files = dir('*.mat');
numFiles = length(files);
for i = 1:numFiles
    fileName = files(i).name;
    a = load(fileName);
    if ~isfield(a, 'var')
        [fileName, ' doesnt have var as a field.']
        var = zeros(1,4);
        var(1) = a.recon.alpha;
        var(2) = a.recon.beta;
        var(3) = a.recon.startPix(1);
        var(4) = a.recon.startPix(2);
        save(fileName,'var', '-v7.3');
        [fileName, ' is over-written.']
    else
        [fileName, ' is fine.']
    end
end

