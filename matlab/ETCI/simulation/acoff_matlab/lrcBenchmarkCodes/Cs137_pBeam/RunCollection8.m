% Run CollectEvents8 on Processed Files (_TRK.mat files) 
LRC = false; % Settings for LRC, if false, Bubbles setting's are used. 
%modified for lawrencium

% To Be located in Data File Folder
% ex1. /global/lr3/acoffer/g4Data/CCD_Tracking/662keV_pBeam/proc_662keV_pBeam

%addpath ./	%not needed
%addpath /global/scratch/bcplimley/multi_angle/MatlabFunctions/
if LRC
    % LRC:
    %addpath ~/ETCI/simulation/acoff_matlab/benchmarkModel/
    % Matlab Scripts Location
    %addpath ~/ETCI/simulation/acoff_matlab/lrcBenchmarkCodes/
    addpath /global/home/users/acoffer/ETCI/simulation/acoff_matlab/benchmarkModel/CCDSpecResponseBenchmark
    workingDir = pwd;
    loadDataPath = workingDir;
	%LRC:
    folderpath = ('/global/lr3/acoffer/g4Data/CCD_Tracking/662keV_pBeam/proc_662keV_pBeam/');  
    %loadDataPath = ('/global/lr3/acoffer/g4Data/CCD_Tracking/662keV_pBeam');
    loadpre = 'Event_pBeam_662KeV_TM_10nmRangCuts_1000nm_LRC_1mil';
    loadsuf = '.mat';
    
    %savepath = '/global/lr3/acoffer/g4Data/CCD_Tracking/662keV_pBeam/proc_662keV_pBeam/';
    %savepre = 'Event_pBeam_662KeV_TM_10nmRangCuts_1000nm_LRC_1mil';
    %savesuf = '.mat';
    %tmp = load('psft.mat');
    % LRC

else
    % Matlab Scripts Location 
    % Bubbles
    addpath /Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoff_matlab/benchmarkModel/CCDSpecResponseBenchmark/
    %addpath /Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoff_matlab/lrcBenchmarkCodes/Cs137_pBeam/
    
    % Bubbles:
    % Folder location 
    folderpath = ('/Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoffer_g4work/CCD_Tracking/temp_data/pBeam_tests_Dec2013/'); 
    loadDataPath = '/Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoffer_g4work/CCD_Tracking/temp_data/pBeam_tests_Dec2013';
    
    loadpre = 'Event_pBeam_662KeV_TM_10nmRangCuts_1000nm_LRC_1mil_';
    loadsuf = '.mat';
    
    %savepath = '/Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoffer_g4work/CCD_Tracking/temp_data/pbeam_tests_dec2013/matlabOutput3/';
    %savepre = '662KeV_pBeam_1000nm_10k_';
    %savesuf = '_TRK.mat';
    %tmp = load('/Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoff_matlab/benchmarkModel/psft_650um.mat');
    %psft = tmp.psft;
end
% Run CollectEvents#.m
E = 0;
out = collectEvents8(loadpre,E, 1);