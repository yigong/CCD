function effCoinc = calculateCoincidentRate(filenamePhrase)

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

effCoinc = 0;
fullECoincEff = 0;
nps = length(fnames)*100000;

for i=1:length(fnames)

    load(fnames{i});
    effCoinc = effCoinc + length(RTrk{1,1});
    
    if 
    
       
end

effCoinc = effCoinc/nps;
