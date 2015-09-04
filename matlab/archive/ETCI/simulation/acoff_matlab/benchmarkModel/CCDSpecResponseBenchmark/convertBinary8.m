function convertBinary8(filenamePhrase, runNPS, depthCordinates)     %phrase characteristic of filename
% convertBinary3('filenamePhrase')

IDphrase = filenamePhrase;
f = dir;
n=1;
%Event = {};

disp('Identifying filenames...')
for i=1:length(dir)
    if ~isempty(strfind(f(i).name,IDphrase))
        fnames{n} = f(i).name;
        s = strfind(fnames{n},'k');
        n=n+1;
    end
end

disp(['Found ',num2str(length(fnames)),' files to be processes'])

% Defaults for Running RunTracks and Geant4TrackHandling
if isempty(runNPS)
    runNPS = 100;
end
if isempty(depthCordinates)
    depthCordinates = [0.002 0.643];
end
    % Load needed matrixs
    %load('../../psft_650um.mat')  % 650um psft, realative path. 
    % Absolute path to file: 
    load('/Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoff_matlab/benchmarkModel/psft_650um.mat')
    %load('../../dEdx_ref.mat')  % Silicon dEdx ref table  % For Hybrid Trk

for i=1:length(fnames)
    disp(['Analyzing',fnames{i},'...'])
    % Start CPU Clock
    t=cputime; %your_operation; cputime-t
    %
    % Convert Binary Track Matrix to Events of Choice
    %
    TrackM = GetCoincidentETrackEvents8(fnames{i}); % Saves events of choice
    %
    % Run Analysis on each Batch of Events
    % 
    [Event] = RunTracks8(TrackM,20000, 'DEdx_ref', psft,depthCordinates); % Saves coincident events
    
    % Save Proccessed data as MAT files
    F = strfind(fnames{i},'.dat');
    sname1 = [fnames{i}(1:(F-1)),'_TRK_processed.mat'];
    sname2 = [fnames{i}(1:(F-1)),'_TRK.mat']
    
    save(sname1,'Event','runNPS');
    disp([' Saving Electron Images and Event Data', sname1, ' to file'])
    
    save(sname2,'TrackM');
    disp([' Saving Electron Track Matrix Data', sname2, ' to file', '\n'])
    
    runTime = (cputime-t)/60;
    disp(['Total time to run ', fnames{i}, ' was ', num2str(runTime), 'min.'])
    
    % Reset Variables
    clear Event TrackM;

    disp(['Done Analyzing ',fnames{i},'...'])

end

