function ARMdist = calculateAllARM(filenamePhrase)
% function ARMdist = calculateAllARM(filenamePhrase)
% phrase characteristic of filename

IDphrase = filenamePhrase;
f = dir;
n=1;

disp('Identifying filenames...')
for i=1:length(dir)
    if ~isempty(strfind(f(i).name,IDphrase))
        fnames{n} = f(i).name;
        s = strfind(fnames{n},'k');
        n=n+1;
    end
end

ARMValues = [];
thetaValues = [];
for i=1:length(fnames)
    disp(['Analyzing',fnames{i},'...'])
    tempHolder = calculateARM(fnames{i}); % Saves calculated events
    tempgARM = tempHolder.gARMdeg;
    
    ARMValues = [ARMValues tempgARM];
    thetaValues = [thetaValues tempHolder.theta];
    disp(['Done Analyzing ',fnames{i},'.../n'])
end


