%run Geant4TrackHandling8 on tracks from 
LRC = true; % Settings for LRC, if false, Bubbles setting's are used. 
%modified for lawrencium

% revised 11/28/12 (and 1/3/13, 5/21/13) after IEEE.
% now lives in /global/scratch/bcplimley/multi_angle/MatlabScripts
% new version: now lives in 

%addpath ./	%not needed
%addpath /global/scratch/bcplimley/multi_angle/MatlabFunctions/
if LRC
    % LRC:
    %addpath ~/ETCI/simulation/acoff_matlab/benchmarkModel/
    % Matlab Scripts Location
    addpath ~/ETCI/simulation/acoff_matlab/lrcBenchmarkCodes/
    addpath ~/ETCI/simulation/acoff_matlab/benchmarkModel/CCDSpecResponseBenchmark/
    %LRC:
    folderpath = ('/global/lr3/acoffer/g4Data/CCD_Tracking/662keV_pBeam/');  
    loadDataPath = ('/global/lr3/acoffer/g4Data/CCD_Tracking/662keV_pBeam/');
    loadpre = 'Mat_662KeV_TM_10nmRangCuts_1000nm_LRC_1mil';
    loadsuf = '.dat';
    
    savepath = '/global/lr3/acoffer/g4Data/CCD_Tracking/662keV_pBeam/proc_662keV_pBeam/';
    savepre = 'Event_pBeam_662KeV_TM_10nmRangCuts_1000nm_LRC_1mil';
    savesuf = '.mat';
    %tmp = load('psft.mat');
    % LRC
    tmp = load('~/ETCI/simulation/acoff_matlab/benchmarkModel/psft_650um.mat');
    psft = tmp.psft;
else
    % Matlab Scripts Location 
    % Bubbles
    addpath /Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoff_matlab/benchmarkModel/CCDSpecResponseBenchmark/
    addpath /Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoff_matlab/lrcBenchmarkCodes/
    
    % Bubbles:
    % Folder location 
    folderpath = ('/Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoffer_g4work/CCD_Tracking/temp_data/pBeam_tests_Dec2013/'); 
    loadDataPath = '/Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoffer_g4work/CCD_Tracking/temp_data/pBeam_tests_Dec2013';
    
    loadpre = 'Mat_662KeV_pBeam_1000nm_10k_';
    loadsuf = '.dat';
    
    savepath = '/Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoffer_g4work/CCD_Tracking/temp_data/pbeam_tests_dec2013/matlabOutput3/';
    savepre = '662KeV_pBeam_1000nm_10k_';
    savesuf = '_TRK.mat';
    tmp = load('/Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoff_matlab/benchmarkModel/psft_650um.mat');
    psft = tmp.psft;
end

placeholdersuf = '_PH.mat';
%which folders to look in? 
% foldernum = {'01','02','03','04','05','06','07','08'};
%foldernum = {'01'};

% Parameter values
%pixelSize = [2.5, 5, 10.5, 20, 30, 40, 50];
pixelSize = 10.5;
%index:       1   2   3     4   5   6   7
% nb = 0.01874;	%base noise level [keV sigma / pixel]
% 1/3/13: only vary noise for 10.5 um.
%noise = [0, 15, 20, 50, 100, 200, 500, 1e3, 2e3];
noise = 18.74;
pixelNoiseMatrix = true;
depthCordinates = [0.002 0.643];
%pixelNoiseMatrix = false(length(pixelSize),length(noise));
%pixelNoiseMatrix(:,1:2) = true; %evaluate all pixel sizes with 0 and 15 eV noise
%pixelNoiseMatrix(3,:) = true;   %evaluate 10.5 um with all noise values

%if we are out of licenses, then return control to shell script
try
    if LRC
        matlabpool local 4;
    else
        % Bubbles
        matlabpool local 1;
    end
catch err
    disp(' ')
    disp('Error:')
    disp(err.identifier)
    disp(err.message)
    disp(' ')
    %pause(2)%*3600)     %[hours]*3600
    if LRC
        exit
    else
        return
    end
    
end

%for j = 1:length(foldernum)     %each folder... (as if we'll ever finish the first one)
    %loadpath = ['/global/scratch/bcplimley/multi_angle/G4batch',foldernum{j},'/'];
    %loadpath = '/global/lr3/acoffer/g4Data/CCD_Tracking/662keV_pBeam/';
    %savepath = '/global/lr3/acoffer/g4Data/CCD_Tracking/662keV_pBeam/proc_662keV_pBeam/';
    
    flist = dir(fullfile(loadDataPath,[loadpre,'*',loadsuf]));
    disp(['Found ',num2str(length(flist)),' files in ',loadDataPath,' at ',datestr(now)])
    
    %files in parallel to avoid conflicts
    parfor i = 1:length(flist)
        try
            %loadDataPath
            %cd loadDataPath;
            G4THworker3(loadDataPath,flist,savepath,savepre,loadpre,loadsuf,savesuf,noise,pixelSize,pixelNoiseMatrix,i,placeholdersuf,psft, depthCordinates);
        catch err
            disp(['Error: ',err.message,' on ',flist(i).name])
        end
    end
%end
matlabpool close;
