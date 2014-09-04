function res = plotDepletionThickness_ESpec(E,plotFlag)
% plotDepletionThickness_ESpec(E,plotFlag)

path(path, '/Users/amycoffer/Research/BeARING/ETI/expBenchmarkModel/depletionCCDSpecsBench')

D = dir;
numFiles= 1;
fileIndex = [];
y = [];
res = {};

% File info: 
% Smallest Depletion to Largest Depletion Thickness energy

%               1   2   3
dThickness = [598 600 650]; % [um]
% Choose an Depletion Thickness to Compare to 
imgSpecNum = 1
 % real nps (check: /mnt/grimdata5/coffer/Spec_benchmark_Sept2012/CCDSpecBatch01 )
 %      /598_DepletionRegion 
 %      /600_Depletion
 %      /650_Depletion 
 %npsCheat = [(38*(5e5)) (29*(5e5)) (30*(5e5)) ];
 % test
 npsCheat = [(2*(5e5)) (1*(5e5)) (1*(5e5)) ];

 
 numel(npsCheat)

% 2 Energies
%fileEnergies = [80.9971  356.017 ];
%gamBranchRatio = [.3406  .6205];  % On 4/12, changed 356keV gamBranchRatio to 0.62 from the incorrect 0.72
%nps = [(79*10e7) (83*(10e7))];

% Set-Up Histograph Parameters
numBins = 200;
maxEnergyBin = 700; % [keV]
binWidth = maxEnergyBin/numBins;
x1 = 0:binWidth:maxEnergyBin;

n = 0;
histTotal = zeros(3,length(x1));
tempEimgTot = [];

for i = 1:length(D);
    
    if( (D(i,1).isdir == 1) && ~strcmp(D(i,1).name , '.') && ~strcmp(D(i,1).name, '..') )
        % Update numFiles and Save fileIndex
        %fileIndex = [fileIndex i]
        %EconTemp = [];
        n = n+1
        
        % Now Have files that are filled with Energies:
        filename = [D(i,1).name '/' 'SummaryEventInfoTRK.mat']
       
        load(filename)
        % 
        % [y(numFiles)] = hist(Etot,x1); % Update Histogram
        
        % Update E Total( CCD Spectrum)
        NTemp = hist(EimgTot, x1);
        
        histTotal(n,:) = NTemp;
        
        if n == imgSpecNum
            tempEimgTot = EimgTot;
        end
    end
end
 
    
% ============   
% === Plot ===
% ============
%plotFlag = 1;
if plotFlag ==1
   
    figure(2)
    % Simulated Cs-137
    %subplot(2,1,2)
    plot(x1, histTotal(1,:)./npsCheat(1),  'b', ...
             x1, histTotal(2,:)./npsCheat(2),  'k', ...
             x1, histTotal(3,:)./npsCheat(3),  'r','LineWidth',2);
    legend('598\mum','600\mum',...
           '650\mum', 'Location','SouthWest');    
    title('Simulated Cs-137 - CCD Spectrum', 'Fontsize',14)
    xlabel('CCD Event Energy [keV]        ','Fontsize',14);
    ylabel('Probability per Incident Photon    ','Fontsize',14);
    
    % Simulation vs Experimental
    % Set-Up Histograph Parameters
    numBins = 200;
    maxEnergyBin = 700; % [keV]
    binWidth = maxEnergyBin/numBins;
    x1 = 0:binWidth:maxEnergyBin;
    y4 = hist(E, x1);


    % Use Number of Events Under the Compton Continuum 
    npsExp = sum(E>100 & E< 700)
    npsSim = sum(tempEimgTot > 100 & tempEimgTot < 670)
    
    
    
    figure(4);
    plot(x1, y4./npsExp, 'b', ...
         x1,histTotal(imgSpecNum,:)./npsSim,'r', ...
         'LineWidth',2);
    legend('Experimental Energy Spectrum', 'Simulated CCD Energy Spectrum');
    title('CCD Total Energy Spectrum','Fontsize',14);  
    xlabel('CCD Event Energy [keV]      ','Fontsize',12);
    ylabel('Event Probability     ','Fontsize',12);
    
    res.normExperimentalSpec = y4./npsExp;
    
    for i=1:n
        
        res.normSimulatedSpec(i,:) = histTotal(i,:)./npsCheat(i);
    end
    
    res.simulatedSpec = histTotal;
    res.npsCheap = npsCheat;
    
    
    %subplot(2,1,1)
    %expHist = hist(E,x1);
    %semilogy(x1,expHist, 'LineWidth',2); 
    %axis([0 x1(length(x1)) 0 5e4]);
   
    title(' CCD Spectrum     ', 'Fontsize',14)
    xlabel('CCD Event Energy [keV]        ','Fontsize',14);
    ylabel('Prob. per incident photon     ','Fontsize',14);
    % Calculate Experimental Peak Ratios
    %for i = 1:length(fileEnergies)
    %    res.exp(i) = PeakFitting(x1,expHist,[(fileEnergies(i)-ROIsigma) (fileEnergies(i)+ROIsigma)]);
    %end 
    %res.expRatio_356_81 = res.exp(1,4).cnts1/res.exp(1,1).cnts1;
    %res.expRatio_81_356 = res.exp(1,1).cnts1/res.exp(1,4).cnts1;

    
    
end

% Save Info Of Interest
sname = ['CCD_SPEC_DepletionLayers.mat'];
save(sname, 'res');

