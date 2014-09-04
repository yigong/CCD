function convertBinary9(filenamePhrase, runNPS, depthCordinates,radDecayFlag, saveFolder)     %phrase characteristic of filename
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
    runNPS = 1e6;
end
if isempty(depthCordinates)
    depthCordinates = [0.002 0.643];
end
if isempty(saveFolder)
    saveFolder = '';
end

if ~exist('radDecayFlag','var')
    % For now, non-Radioactive Decay Sources are default
    radDecayFlag = false;
end

    % Load needed matrixs
    %load('../../psft_650um.mat')  % 650um psft, realative path. 
    %load('../../dEdx_ref.mat')     % dEdx ref, relative path. 
    % Absolute path to file: 
    load('/Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoff_matlab/benchmarkModel/psft_650um.mat')
    load('/Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoff_matlab/benchmarkModel/dEdx_ref.mat')  % Silicon dEdx ref table  % For Hybrid Trk

for i=1:length(fnames)
    disp(['Analyzing',fnames{i},'...'])
    % Start CPU Clock
    t=cputime; %your_operation; cputime-t
    %
    % Convert Binary Track Matrix to Events of Choice
    %
    TrackM = GetCoincidentETrackEvents9(fnames{i},true); % Saves events of choice (filename, radDecayFlag)
    disp(['Analyzing ', num2str(length(TrackM)), ' possible electron tracks '])
    %
    % Can Load Saves Segmented Tracks at this point too..
    %load(fnames{i}); % Load a data structure with one TrackM{1,i} per tracked Event
    
    % Run Analysis on each Batch of Events
    
    [Event] = RunTracks9(TrackM,20000, 'DEdx_ref', psft,depthCordinates,radDecayFlag); % Saves Proccessed events
    
    % Save Proccessed data as MAT files
    F = strfind(fnames{i},'.dat');
    sname1 = [saveFolder,'/',fnames{i}(1:(F-1)),'_TRK_processed.mat'];
    sname2 = [saveFolder,'/',fnames{i}(1:(F-1)),'_TRK.mat'];
    
    save(sname1,'Event','runNPS');
    disp([' Saving Electron Images and Event Data ', sname1, ' to file \n'])
    
    save(sname2,'TrackM');
    disp([' Saving Electron Track Matrix Data ', sname2, ' to file \n'])
    
    runTime = (cputime-t)/60;
    disp(['Total time to run ', fnames{i}, ' was ', num2str(runTime), 'min. \n'])
    
    % Reset Variables
    clear Event TrackM;

    disp(['Done Analyzing ',fnames{i},'...'])

end

