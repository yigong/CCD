%run Geant4TrackHandling7 on tracks from /mnt/grimdata5/coffer/eT_MultiAngle_Si_Sept2012/Batch01/

%modified for lawrencium

% revised 11/28/12 (and 1/3/13, 5/21/13) after IEEE.
% now lives in /global/scratch/bcplimley/multi_angle/MatlabScripts

addpath ./	%not needed
addpath /global/scratch/bcplimley/multi_angle/MatlabFunctions/

%which folders to look in?
% foldernum = {'01','02','03','04','05','06','07','08'};
foldernum = {'01'};

%what parameter values?
pixelSize = [2.5, 5, 10.5, 20, 30, 40, 50];
%index:       1   2   3     4   5   6   7
% nb = 0.01874;	%base noise level [keV sigma / pixel]

% 1/3/13: only vary noise for 10.5 um.
noise = [0, 15, 20, 50, 100, 200, 500, 1e3, 2e3];
pixelNoiseMatrix = false(length(pixelSize),length(noise));
pixelNoiseMatrix(:,1:2) = true; %evaluate all pixel sizes with 0 and 15 eV noise
pixelNoiseMatrix(3,:) = true;   %evaluate 10.5 um with all noise values


%5/21/13: I want to restore full matrix. But run algorithm on existing segmented tracks first...
%see logbook 11/28 for noise scaling.
%{
noise = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1e3, 2e3, 5e3, 1e4, 2e4];
%index:  1  2  3  4   5   6   7   8    9   10   11   12   13   14   15

% see logbook 11/28. smaller pixels scale with area; larger pixels scale with length.
pixelNoiseMatrix = true(length(pixelSize),length(noise));
%right blacked-out part
pixelNoiseMatrix(1, end-6:end)  = false;
pixelNoiseMatrix(2, end-4:end)  = false;
pixelNoiseMatrix(3, end-2:end)  = false;
pixelNoiseMatrix(4, end)        = false;
pixelNoiseMatrix(5, end)        = false;
%left blacked-out part
pixelNoiseMatrix(2, 2)          = false;
pixelNoiseMatrix(3, 2:4)        = false;
pixelNoiseMatrix(4, 2:6)        = false;
pixelNoiseMatrix(5, 2:6)        = false;
pixelNoiseMatrix(6, 2:7)        = false;
pixelNoiseMatrix(7, 2:7)        = false;
%}

%GRIM:
%loadpath = '/mnt/grimdata5/coffer/eT_MultiAngle_Si_Sept2012/eT_MultiAngle_Si_Batch01/';
%LRC:
%loadpath moved inside "for j" loop
loadpre = 'Mat_MultiAngle_CCD_eTracks_662_500k_';
loadsuf = '_TRK.mat';

%GRIM:
%savepath = '/mnt/grimdata5/plimley/multi_angle_2012/Batch01/';
%LRC:
%savepath moved inside "for j" loop
savepre = 'MultiAngle_DT_';
savesuf = '.mat';
placeholdersuf = '_PH.mat';

tmp = load('psft.mat');
psft = tmp.psft;

%if we are out of licenses, then return control to shell script
try
    matlabpool local 4;
catch err
    disp(' ')
    disp('Error:')
    disp(err.identifier)
    disp(err.message)
    disp(' ')
    pause(2*3600)     %[hours]*3600
    exit
end

for j = 1:length(foldernum)     %each folder... (as if we'll ever finish the first one)
    loadpath = ['/global/scratch/bcplimley/multi_angle/G4batch',foldernum{j},'/'];
    savepath = ['/global/scratch/bcplimley/multi_angle/DTbatch',foldernum{j},'/'];
    
    flist = dir(fullfile(loadpath,[loadpre,'*',loadsuf]));
    disp(['Found ',num2str(length(flist)),' files in ',loadpath,' at ',datestr(now)])
    
    %files in parallel to avoid conflicts
    parfor i = 1:length(flist)
        try
            G4THworker(loadpath,flist,savepath,savepre,loadpre,loadsuf,savesuf,noise,pixelSize,pixelNoiseMatrix,i,placeholdersuf,psft);
        catch err
            disp(['Error: ',err.message,' on ',flist(i).name])
        end
    end
end
matlabpool close;
