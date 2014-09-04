% Plotting 
% New Tests of CCD-Spec
% 
%function collectEvents1(filenamePhrase,E,plotFlag)     
%path(path,'/Users/amycoffer/Research/BeARING/ETI/expBenchmarkModel/CCDSpecResponseBenchmark/CCDSpec_September2012')

% Load Collected Events Output file
%if ~exist('EimgTot','var')
%    load('CCDSpecBench_Experiment_and_Simulation.mat')
%end

plotFlag = 1;
if ~exist('E','var')
    addpath /Users/amycoffer/Dropbox/Uni_Research/ETI/ETCI/simulation/acoffer_g4work/data/experimentalData
    % Beam Source: (verfiy?)
    load('coincBench_CCD_ESpec_trackE20_0.mat');
    % pt. source
    %load('xyzlist_cs137_q1 (1).mat'); % pt source E data
    %E = 0;
end

% Set-Up Histograph Parameters
numBins = 350;
maxEnergyBin = 700; % [keV]
binWidth = maxEnergyBin/numBins;
y1 = 0:binWidth:maxEnergyBin;

% Normaized Energy Spec of Simulation and Experiment to the total Number of
% Counts above 100 keV
Emin = 100;     % [keV]
Emax = 700;

% Use Number of Events Under the Compton Continuum 
npsExp = sum(E>Emin & E< Emax)
npsSim1 = sum(EimgTot > Emin & EimgTot < Emax)

% for i=1:length(E)
%     if E(i) > Emin
%         Eprime = [Eprime E(i)];
%         npsExp = npsExp+1;
%     end
% end

% Create Modified EimgTot by adding some rate of track overlap
modEimgTot = EimgTot(isfinite(EimgTot(1,:)));
overlap = 0.015; 
numDoubles = round(overlap*npsSim1);
npsSim = npsSim1 + numDoubles;
R = randi(length(modEimgTot), [2,numDoubles]);
tempEimgtot = modEimgTot(R(1,:))+modEimgTot(R(2,:));
EimgTot2 = cat(2,modEimgTot,tempEimgtot); 
npsSim = sum(EimgTot2 > Emin & EimgTot2 < Emax);

[x0,y1] = hist(EimgTot2,y1);
[x1,y0] = hist(Edep,y1);
[x2,y2] = hist(EimgContained, y1);
[x3,y3] = hist(EimgNonContained, y1);
[x4,y4] = hist(E, y1);

% Save Histogram Info

%sname = ['SummaryInfo', filenamePhrase];
%save(sname, 'y0', 'y1', 'y2', 'y3', 'y4', 'x1', 'x2', 'x3', 'x4');

%sname2 = ['SummaryEventInfo', filenamePhrase];
%save(sname2, 'Etot', 'Edep', 'Eesc', 'Econtained','Enoncontained', 'multi');


% Plot Histogram Info
if plotFlag == 1
    plot(y1, x0, 'g', y1,x2,'b', y1,x3,'r')
    
    % Calculate Initial Number of Gamma's
    %nps = 101*500000 %4*10*10^6+2*5*10^6; %Total Number of Starting Photons
    %nps2 = length(E);
    nps = 20*(1e6);
    figure(1)
    % Simulation
    %subplot(2,1,2);
    plot(y1, x0./nps, 'k', ...
         y1, x2./nps, 'b', ...
         y1, x3./nps, 'r','LineWidth',2);

    legend('Total Energy Spectrum','Cointained Electrons','Non-contained Electrons');
    title('Simulated CCD Energy Spectrum', 'Fontsize',14)
    xlabel('CCD Event Energy [keV]     ','Fontsize',12);
    ylabel('Event Probability     ','Fontsize',12);
    %axes('Fontsize',20)
    
    % Experiment and Simulation Total CCD-Energy
    %subplot(2,1,1);
    figure(2)
%    plot(y1,x4./npsExp, 'b', ...
%         y1,x0./npsSim,'r', ...
    plot(y1,x4, 'b', ...
         y1,x0,'r', ... 
            'LineWidth',2);
    legend('Experimental Energy Spectrum', 'Simulated CCD Energy Spectrum');
    title('CCD Total Energy Spectrum','Fontsize',14);  
    xlabel('CCD Event Energy [keV]      ','Fontsize',12);
    ylabel('Event Probability     ','Fontsize',12);

end