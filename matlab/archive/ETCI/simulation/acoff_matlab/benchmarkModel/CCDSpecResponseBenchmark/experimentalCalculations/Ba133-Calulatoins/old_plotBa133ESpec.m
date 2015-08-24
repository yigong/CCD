function res = plotBa133ESpec(E,plotFlag)
% plotBa133ESpec(E,plotFlag)
% Proper Path:
% /Users/amycoffer/Research/BeARING/ETI/expBenchmarkModel/CCDSpecResponseBenchmark
path(path, '/Users/amycoffer/Research/BeARING/ETI/expBenchmarkModel/newBa133_Simulations');

D = dir;
numFiles= 5;
ROIsigma = 6;
fileIndex = [];
y = [];
res = {};

% File info: 
% Low energy to high energy

% 5 Energies      [   1         2         3         4         5     ];
 fileEnergies =   [80.9971  276.398   302.853   356.017   383.851];
 gamBranchRatio = [  .3406     .07164    .1833     .6205     .0894];
 % old: npsCheat = [(79*(10e7)) (12*(10e7)) (12*(10e7)) (83*(10e7)) (13*(10e7))]
 % new nsp:
 %%npsCheat = [(79*(10e7)) (53*(10e7)) (51*(10e7)) (83*(10e7)) (46*(10e7))]
 %npsCheat = [(*(10e7)) (*(10e7)) (*(10e7)) (*(10e7)) (*(10e7))]
 npsCheat = [5.7e+09 3.1e+09 4.9e+09 4.9e+09 17500000]; 
 dirOrder = [ 2 3 4 5 1];
 numel(npsCheat)

% 2 Energies
%fileEnergies = [80.9971  356.017 ];
%gamBranchRatio = [.3406  .6205];  % On 4/12, changed 356keV gamBranchRatio to 0.62 from the incorrect 0.72
%nps = [(79*10e7) (83*(10e7))];

% Set-Up Histograph Parameters
numBins = 200;
maxEnergyBin = 400; % [keV]
binWidth = maxEnergyBin/numBins;
x1 = 0:binWidth:maxEnergyBin;

% Temp Holders
histTotal = zeros(1,length(x1));
histEIC = zeros(1,length(x1));
histEINC = zeros(1,length(x1));

n = 0;

for i = 1:length(D);
    % Build Ba-133 Spectrum 
    if( (D(i,1).isdir == 1) && ~strcmp(D(i,1).name , '.') && ~strcmp(D(i,1).name, '..') )
        % Update numFiles and Save fileIndex
        %fileIndex = [fileIndex i]
        EconTemp = [];
        n = n+1
        
        % Now Have files that are filled with Energies:
        filename = [D(i,1).name '/' 'SummaryEventInfoTRK.mat']
        load(filename)
        % 
        % [y(numFiles)] = hist(Etot,x1); % Update Histogram
        
        % Update E Total( CCD Spectrum)
        histTemp = hist(EimgTot, x1);
        normHistTemp = histTemp*gamBranchRatio(dirOrder(n))/npsCheat(dirOrder(n));
        histTotal = histTotal + normHistTemp;
        
        % Update EimgNonContained
        histTempEIC = hist(EimgContained, x1);
        normHistTempEIC = histTempEIC*gamBranchRatio(dirOrder(n))/npsCheat(dirOrder(n));
        histEIC = histEIC + normHistTempEIC;
        % Update EimgContained
        histTempEINC = hist(EimgNonContained, x1);
        normHistTempEINC = histTempEINC*gamBranchRatio(dirOrder(n))/npsCheat(dirOrder(n));
        histEINC = histEINC +normHistTempEINC;
        
        % Clear TempHists
        histTemp = []; histTempEIC = []; histTempEINC = [];
        
    end
end

% Peak Fitting and Peak Ratio Calculations:
%plot(x1,histTotal)
%res1 = PeakFitting(b,n,81)
for i=1:length(fileEnergies)
    res.sim(i) = PeakFitting(x1,histTotal,[(fileEnergies(i)-ROIsigma) (fileEnergies(i)+ROIsigma)]);
    %res2 = PeakFitting(b,n,356)
    %res2 = PeakFitting(x1,histTotal,[350 361])
end


    
% Hist of interest
    res.histTotal = [histTotal x1];
    
% Calculate Experimental Peak Ratios
    expHist = hist(E,x1);
    % Calculate Peak Info For all 5 Energies
    for i = 1:length(fileEnergies)
        res.exp(i) = PeakFitting(x1,expHist,[(fileEnergies(i)-ROIsigma) (fileEnergies(i)+ROIsigma)]);
    end 
    
    % Calculate Simulation Ratios
    % Calculate Simulated Peak Ratio:
%for i = 1:length(fileEnergies)
    %                   1       2       3       4       5
    % fileEnergies = [80.9971 276.398 302.853 356.017 383.851];

    res.simRatio_356_81 = (res.sim(1,4).cnts1)/(res.sim(1,1).cnts1);
    res.simRatio_356_276 = res.sim(1,4).cnts1/res.sim(1,2).cnts1;
    res.simRatio_356_303 = res.sim(1,4).cnts1/res.sim(1,3).cnts1;
    res.simRatio_356_384 = res.sim(1,4).cnts1/res.sim(1,5).cnts1;

    res.simRatio_81_276 = res.sim(1,1).cnts1/res.sim(1,2).cnts1;
    res.simRatio_81_303 = res.sim(1,1).cnts1/res.sim(1,3).cnts1;
    res.simRatio_81_356 = res.sim(1,1).cnts1/res.sim(1,4).cnts1;
    res.simRatio_81_384 = res.sim(1,1).cnts1/res.sim(1,5).cnts1;
    
    % Calculate Experimental Ratios
    res.expRatio_356_81 = res.exp(1,4).cnts1/res.exp(1,1).cnts1;
    res.expRatio_356_276 = res.exp(1,4).cnts1/res.exp(1,2).cnts1;
    res.expRatio_356_303 = res.exp(1,4).cnts1/res.exp(1,3).cnts1;
    res.expRatio_356_384 = res.exp(1,4).cnts1/res.exp(1,5).cnts1;
    
    res.expRatio_81_276 = res.exp(1,1).cnts1/res.exp(1,2).cnts1;
    res.expRatio_81_303 = res.exp(1,1).cnts1/res.exp(1,3).cnts1;
    res.expRatio_81_356 = res.exp(1,1).cnts1/res.exp(1,4).cnts1;
    res.expRatio_81_384 = res.exp(1,1).cnts1/res.exp(1,5).cnts1;
    
    % Calculate Percent-Difference 
    % Between Experimental and Simulation Values of Peak Ratios
    
    % === 356 ===
    % 356/81
    res.percentDiff_356_81 = ((res.simRatio_356_81)-(res.expRatio_356_81))/res.expRatio_356_81;
    % 356/276
    res.percentDiff_356_276 = ((res.simRatio_356_276)-(res.expRatio_356_276))/res.expRatio_356_276;
    % 356/303
    res.percentDiff_356_303 = ((res.simRatio_356_303)-(res.expRatio_356_303))/res.expRatio_356_303;
    % 356/384
    res.percentDiff_356_384 = ((res.simRatio_356_384)-(res.expRatio_356_384))/res.expRatio_356_384;
    
    
    % === 81 ===
    % 81 /276
    res.percentDiff_81_276 = ((res.simRatio_81_276)-(res.expRatio_81_276))/res.expRatio_81_276;
    % 81/303
    res.percentDiff_81_303 = ((res.simRatio_81_303)-(res.expRatio_81_303))/res.expRatio_81_303;
    % 81/356
    res.percentDiff_81_356 = ((res.simRatio_81_356)-(res.expRatio_81_356))/res.expRatio_81_356;
    % 81/384
    res.percentDiff_81_384 = ((res.simRatio_81_384)-(res.expRatio_81_384))/res.expRatio_81_384;
    
    
    
% ============   
% === Plot ===
% ============
%plotFlag = 1;
if plotFlag ==1
   
    figure(4)
    % Experimental Ba-133
    subplot(2,1,2)
    plot(x1, histEIC,  'b', ...
             x1, histTotal,  'k','LineWidth',2); % ...
             %x1, histEINC,  'r', %semilogy(x1, histTotal)
    legend('Total Energy Spectrum','Cointained Electrons',...
           'Location','SouthWest');    
            % 'Non-contained Electrons'
    title('Simulated Ba-133 - CCD Spectrum', 'Fontsize',14)
    xlabel('CCD Event Energy [keV]','Fontsize',14);
    ylabel('Event Probability     ','Fontsize',14);
    % Simulation Ba-133
    subplot(2,1,1)
    expHist = hist(E,x1);
    plot(x1,expHist, 'LineWidth',2); 
    axis([0 x1(length(x1)) 0 5e4]);
   
    title('Experimental Ba-133 - CCD Spectrum', 'Fontsize',14)
    xlabel('CCD Event Energy [keV]','Fontsize',14);
    ylabel('Counts    ','Fontsize',14);
    % Calculate Experimental Peak Ratios
    %for i = 1:length(fileEnergies)
    %    res.exp(i) = PeakFitting(x1,expHist,[(fileEnergies(i)-ROIsigma) (fileEnergies(i)+ROIsigma)]);
    %end 
    %res.expRatio_356_81 = res.exp(1,4).cnts1/res.exp(1,1).cnts1;
    %res.expRatio_81_356 = res.exp(1,1).cnts1/res.exp(1,4).cnts1;

    
    
end

% Save Info Of Interest
sname = ['Ba133-HistogramData.mat'];
save(sname, 'res');



% % Set-Up Histograph Parameters
% numBins = 200;
% maxEnergyBin = 700; % [keV]
% binWidth = maxEnergyBin/numBins;
% y1 = 0:binWidth:maxEnergyBin;
% 
% %[x0,y1] = hist(Etot,200);
% [x1,y0] = hist(Edep,y1);
% [x2,y2] = hist(Econtained, y1);
% [x3,y3] = hist(Enoncontained, y1);
% [x4,y4] = hist(E, y1);
%     
% % Calculate Initial Number of Gamma's
%     nps = 101*500000; %4*10*10^6+2*5*10^6; %Total Number of Starting
%     Photons
%  
%     
%     figure(2)
%     subplot(2,1,2);plot(y1, x1./nps, 'k', ...
%          y1, x2./nps, 'b', ...
%          y1, x3./nps, 'r','LineWidth',2);
% 
%     legend('Total Energy Spectrum','Cointained Electrons','Non-contained Electrons');
%     title('Simulated CCD Spectrum of Experiment', 'Fontsize',14)
%     xlabel('CCD Event Energy [keV]','Fontsize',14);
%     ylabel('Event Probability','Fontsize',14);
%     %axes('Fontsize',20)
%     nps2 = length(E);
%     subplot(2,1,1);plot(y1,x4./nps2,'b','LineWidth',2);
%     title('Experimental CCD Total Energy Spectrum','Fontsize',14);  
%     xlabel('CCD Event Energy [keV]','Fontsize',14);
%     ylabel('Event Probability [arb units]','Fontsize',14);