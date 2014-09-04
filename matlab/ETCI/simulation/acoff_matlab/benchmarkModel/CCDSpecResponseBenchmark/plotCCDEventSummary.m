% Plot 
%function collectEvents1(filenamePhrase,E,plotFlag)     
path(path,'/Users/amycoffer/Research/BeARING/ETI/expBenchmarkModel/CCDSpecResponseBenchmark/finalSimulationBeamVectorSpec')

if ~exist('EimgTot','var')
    load('CCDSpecBench_Experiment_and_Simulation.mat')
end

plotFlag = 1;
if ~exist('E','var')
    E = 0;
end

% Set-Up Histograph Parameters
numBins = 440;
maxEnergyBin = 680; % [keV]
binWidth = maxEnergyBin/numBins;
y1 = 0:binWidth:maxEnergyBin;

% Normaized Energy Spec of Simulation and Experiment to the total Number of
% Counts above 100 keV
Emin = 100;     % [keV]

% Use Number of Events Under the Compton Continuum 
npsExp = sum(E>100 & E< 700)
npsSim = sum(EimgTot > 100 & EimgTot < 700)

% for i=1:length(E)
%     if E(i) > Emin
%         Eprime = [Eprime E(i)];
%         npsExp = npsExp+1;
%     end
% end
       
[x0,y1] = hist(EimgTot,y1);
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
    nps = 101*500000 %4*10*10^6+2*5*10^6; %Total Number of Starting Photons
    %nps2 = length(E);
    
    figure(2)
    % Simulation
    subplot(2,1,2);plot(y1, x0./nps, 'k', ...
         y1, x2./nps, 'b', ...
         y1, x3./nps, 'r','LineWidth',2);

    legend('Total Energy Spectrum','Cointained Electrons','Non-contained Electrons');
    title('Simulated CCD Energy Spectrum', 'Fontsize',14)
    xlabel('CCD Event Energy [keV]     ','Fontsize',12);
    ylabel('Event Probability     ','Fontsize',12);
    %axes('Fontsize',20)
    
    % Experiment and Simulation Total CCD-Energy
    subplot(2,1,1);
    plot(y1,x4./npsExp, 'b', ...
         y1,x0./npsSim,'r', ...
         'LineWidth',2);
    legend('Experimental Energy Spectrum', 'Simulated CCD Energy Spectrum');
    title('CCD Total Energy Spectrum','Fontsize',14);  
    xlabel('CCD Event Energy [keV]      ','Fontsize',12);
    ylabel('Event Probability     ','Fontsize',12);

end