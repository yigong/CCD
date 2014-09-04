% Experimental Estimates for 
% Coincident Experimental Setup:
%       1.) Photopeak to Compton Ratio (Knoll (3rdEd) Def. page 428)
%       2.) Interation Rate estimate
%       3.) Coincident Rate estimate

% cross section data
addpath /Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoff_matlab/benchmarkModel/CCDSpecResponseBenchmark/experimentalCalculations

load('Si_XCOM_10keV-10MeV')
% Experimental Data
if ~exist('E','var')
    addpath /Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoffer_g4work/data/experimentalData
    % Coincident Experimental Setup CCD Energy Energy
    load('coincBench_CCD_ESpec_trackE20_0.mat')  % Loads just energy vector E.
end

%xcom_PE662 = interp1(PhotonEnergy, PhotoelAbsorb,0.662)
xcom_PE662_pchip = interp1(PhotonEnergy1, PhotoElAbsorb,0.662,'pchip') % [cm^2/g]
xcom_Total662_pchip = interp1(PhotonEnergy1, TotWCoherent,0.662,'pchip')

rho_Si = 2.330% [g/cm^3] 
% Other Experimental Values:
sourceStrength =  7.216e7; % [photons/sec]
Al_Transmission = 0.9324;
BeamSpot = 0.444;
% CCD Specs
CCD_Thickness = .650/10; % [cm]
CCD_relativeThickness = 0.174; % [cm], 650um CCD with 22 degree from perp of CCD face. 
CCD_rhoXThickness = rho_Si*CCD_relativeThickness;

P_PE662keV = xcom_PE662_pchip*CCD_rhoXThickness;
P_Total662keV = xcom_Total662_pchip*CCD_rhoXThickness;

PhotoEAbsorb_to_Total =  P_PE662keV/P_Total662keV

numPer500k = PhotoEAbsorb_to_Total *(5e5);
numPer500k_50perEcontained = numPer500k*.5;
numPer500k_20perEcontained = numPer500k*.2;
numPer500k_10perEcontained = numPer500k*.10;

% Experimental Comparison
% Normaized Energy Spec of Simulation and Experiment to the total Number of
% Counts above 100 keV
Emin = 100;     % [keV]
Emax = 700;

% Use Number of Events Under the Compton Continuum 
npsExp = sum(E>Emin & E< Emax)
%npsSim1 = sum(EimgTot > Emin & EimgTot < Emax)

% Histogram Experiment
% Set-Up Histograph Parameters
numBins = 350;
maxEnergyBin = 700; % [keV]
binWidth = maxEnergyBin/numBins;
x1 = 0:binWidth:maxEnergyBin;

% Peak Fitting REI
roiMin = 640; % [keV]
roiMax = 680;

% Calculate Experimental Peak Ratios
    expHist = hist(E,x1);

% Calculate Peak Info for 662keV Peak
    res.exp = PeakFitting(x1,expHist,[roiMin  roiMax]);

    % Simulation
    if exist('EimgTot','var')
        simHist = hist(EimgTot,x1);
        res.sim = PeakFitting(x1,simHist,[roiMin  roiMax]);
    end