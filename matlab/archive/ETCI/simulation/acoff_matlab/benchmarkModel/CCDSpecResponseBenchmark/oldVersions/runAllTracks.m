function runAllTracks(filenamePhrase)     %phrase characteristic of filename

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

for i=1:length(fnames)
    disp(['Analyzing',fnames{i},'...'])
    RunTracks(fnames{i},20000); % Saves coincident events
    disp(['Done Analyzing ',fnames{i},'...'])

end

