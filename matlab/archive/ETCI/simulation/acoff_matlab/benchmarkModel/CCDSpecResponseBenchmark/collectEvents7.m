function out = collectEvents7(filenamePhrase,E, plotFlag)
% out = collectEvents6(filenamePhrase,E, plotFlag)
% [Etot, Edep, EescOld, EescNew, Econtained,Enoncontained,
% multi,error1,error2] = collectEvents2(filenamePhrase,E, plotFlag)
% [Etot, Edep, Eesc, Econtained,Enoncontained, multi] =
% collectEvents2(filenamePhrase,E) 
% filename == _TRK.mat file with Event{} Structure.%phrase characteristic
% of filename
% E = CCD Benchmarking Experiment Tracks Energy Matrix
% plotFlag: 1 if you want histograms plotted.
% saveflag: currently saves two summary files: 
%   SummaryInfo_filenamePhrase - Binned Energy Info
%   SummaryEventInfo_filenamePhase - Energy Info non-binned
% Current as of: 3/14/12 ^^ /\ /\
% -----------------------
% \/ \/
% 10/10/2012:
% collectEvents7.m [mod] (removed ARM and Ge-Det collection loops)
% Modified Version for depletion CCD-Spec Analysis: 
% Only Collects the CCD Spectral Info
                                            
IDphrase = filenamePhrase;
f = dir;
n=1;%
Event = {};
out = {};
primaryPhotons = 500000; % Number of primary photons for each file
badCCD = 0; 

disp('Identifying filenames...')
for i=1:length(dir)
    if ~isempty(strfind(f(i).name,IDphrase))
        fnames{n} = f(i).name;
        s = strfind(fnames{n},'k');
        n=n+1;
    end
end

%ARMValues = [];
%thetaValues = [];
Etot = [];
EimgTot = [];
Edep = [];
EescOld = [];
EescNew = [];
multi = [];
EimgNonContained = [];
EimgContained = [];

% Some Error Tracking
%    error0 = [];
%    error1 = [];
    error2 = [];
    errorCollection = [];   

% Index counter for Number of sucsseful Event files opened...
k= 0;
nps = 0;

for i=1:length(fnames)
    tic
    disp(['Analyzing',fnames{i},'...'])

    try
        % Load i Files of Events
        load(fnames{i});
        %[tempEtot,tempEimgTot, tempEdep,tempEescOld, tempEescNew, tempMulti, tempError1,tempError2] = collectCCDEnergies2(fnames{i}); % Saves calculated events
        if exist('Event','var')
            res = collectCCDEnergies7(Event); % Saves calculated events
            k = k+1;
        else
            disp(['Skipping Analysis of ', fnames{i}, '... Does not contain proper Event Structure... ']);
            continue
            toc
        end
        if exist('runNPS','var')
            nps = nps + runNPS;
        else
            disp(['Event number' num2str(i) ' does not contain runNPS value, will calulate nps bases on primaryPhotons variable which is currently set to ' num2str(primaryPhotons)])
        end            
            
    catch e
        disp(['Skipping Analysis of ', fnames{i}, '... Something errored out..  may not contain proper Event Structure... ']);
        error0 = [error0 e];
        continue
        toc
    end
    

    
    try      
        tic
       
        Etot = [Etot res.Etot];
        EimgTot = [EimgTot res.EimgTot];
        Edep = [Edep res.Edep];
        
        EescOld = [EescOld res.EescOld];
        EescNew= [EescNew res.EescNew];
        multi = [multi res.multi];
        
        EimgNonContained = [EimgNonContained res.EimgNonContained];
        EimgContained = [EimgContained res.EimgContained];
       
        %error1 = [error1 res.error1];
        error2 = [error2 res.error2];
        
        disp(['Done Analyzing ',fnames{i},'...'])
        clear res;
        toc
        
    catch e 
        disp(['Skipping Analysis of ', fnames{i}, '... Bad CCD Energy Information... ']);
        
        badCCD = badCCD+1;
        errorCollection = [errorCollection e];
        toc
        continue
        
    end
 
    clear res 
    clear Event
    
end

% Calculate Number of Initial Photons
if nps<1
    nps = primaryPhotons.*k;
end


% Print Interesting Values:
disp(['Length of geant4 Etot ' num2str(length(Etot)) ' while there are ' num2str(sum(~isnan(Etot))) ' are real numbers.'])
disp(['Length of EimgTot ' num2str(length(EimgTot)) ' while there are ' num2str(sum(~isnan(EimgTot))) ' are real numbers.'])
disp(['There are a recorded ' num2str(badCCD) ' times that a file failed to be collected properly.' ])

if plotFlag ==1
    % Plot Energy Spectrum of Interest
    disp(['Plotting ', filenamePhrase, '...'])
    tic
end



% Set-Up Histograph Parameters
numBins = 300;
maxEnergyBin = 680; % [keV]
binWidth = maxEnergyBin/numBins;
y1 = 0:binWidth:maxEnergyBin;
[x1,y0] = hist(EimgTot,y1);
[x2,y2] = hist(EimgContained, y1);
[x3,y3] = hist(EimgNonContained, y1);
[x4,y4] = hist(E, y1);

% Save Histogram Info
sname = ['SummaryInfo', filenamePhrase];
save(sname, 'y0', 'y1', 'y2', 'y3', 'y4', 'x1', 'x2', 'x3', 'x4');

% Save Histogramable Infos
sname2 = ['SummaryEventInfo', filenamePhrase];
save(sname2, 'Etot', 'EimgTot','Edep', 'EescOld','EescNew', 'EimgContained','EimgNonContained', ...
            'multi', 'nps');

% CCD Info Out
out.Etot = Etot;
out.EimgTot = EimgTot;
out.Edep = Edep;
out.EescOld = EescOld;
out.EescNew = EescNew;
out.EimgContained = EimgContained;
out.EimgNonContained = EimgNonContained;
out.multi = multi;
out.k = k; % num files run
out.nps = nps;
% Errors
%out.errorType0 = error0;
%out.errorType1 = error1;
out.errorType2 = error2;
out.errorCollection = errorCollection;
if exist('runNPS', 'var')
    out.runNPS = runNPS;
end
% Save Histogram Info
sname = ['SummaryInfo', filenamePhrase];
save(sname, 'out','y0', 'y1', 'y2', 'y3', 'y4', 'x1', 'x2', 'x3', 'x4');


% Plot Histogram Info
if plotFlag == 1
    %plot(y1, x1, 'g', y1,x2,'b', y1,x3,'r')
    
    % Calculate Initial Number of Gamma's
    %nps = 1; %100*10*10^6 %4*10*10^6+2*5*10^6; %Total Number of Starting Photons
    tic
    
    figure(3)
    subplot(2,1,2);plot(y1, x1./nps, 'k', ...
         y1, x2./nps, 'b', ...
         y1, x3./nps, 'r','LineWidth',2);

    legend('Total Energy Spectrum','Cointained Electrons','Non-contained Electrons', 'Location','NorthWest');
    title('Simulated CCD Energy Spectrum of Benchmarking Experiment', 'Fontsize',14)
    xlabel('CCD Event Energy [keV]','Fontsize',12);
    ylabel('Coincident Event Probability','Fontsize',12);
    %axes('Fontsize',20)
    nps2 = length(E);
    subplot(2,1,1);
    plot(y1,x4./nps2,'b', ...
         y1, x1./nps2, 'r', 'LineWidth',2);
         %y1, x1./length(Edep), 'r', 'LineWidth',2); %Normalize Ba-133 Spec
         % with total number of events in spectrum. 
    legend('Experimental Energy Spectrum', 'Simulated CCD Energy Spectrum', 'Location','NorthWest');
    title('Experimental CCD Total Energy Spectrum','Fontsize',14);  
    xlabel('CCD Event Energy [keV]','Fontsize',12);
    ylabel('Event Probability [arb units]','Fontsize',12);    
    toc
end

