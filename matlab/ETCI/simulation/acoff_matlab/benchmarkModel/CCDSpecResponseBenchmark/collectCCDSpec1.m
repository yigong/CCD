function [Etot, Edep, Eesc, Econtained,Enoncontained, multi] = collectCCDSpec1(filenamePhrase)     
% [Etot, Edep, Eesc, Econtained,Enoncontained, multi] =
% collectEvents1(filenamePhrase,E) 
% filename == _TRK.mat file with Event{} Structure.%phrase characteristic of filename
% E = CCD Benchmarking Experiment Tracks Energy Matrix
% plotFlag: 1 if you want histograms plotted.
% saveflag: currently saves two summary files: 
%   SummaryInfo_filenamePhrase - Binned Energy Info
%   SummaryEventInfo_filenamePhase - Energy Info non-binned
                                            
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

%ARMValues = [];
%thetaValues = [];
Etot = [];
Edep = [];
Eesc = [];
multi = [];

for i=1:length(fnames)
    disp(['Analyzing',fnames{i},'...'])
    
    
    [tempEtot,tempEdep,tempEesc, tempMulti] = collectCCDEnergies5(fnames{i}); % Saves calculated events
    
    
    Etot = [Etot tempEtot];
    Edep = [Edep tempEdep];
    Eesc = [Eesc tempEesc];
    multi =[multi tempMulti];
    
    tempEtot= []; tempEdep=[]; tempEesc=[]; tempMulti=[];
    
    disp(['Done Analyzing ',fnames{i},'...'])
end

if plotFlag ==1
    % Plot Energy Spectrum of Interest
    disp(['Plotting ', filenamePhrase, '...'])
end

% Determine Energy Spec. of Cointained and Non-Cointained Electrons
Econtained = [];
Enoncontained = [];

 for i=1:length(Etot)
     if Eesc(i) >0
         Enoncontained = [Enoncontained Edep(i)];
     else
         Econtained = [Econtained Edep(i)];
     end
 end
 
% Set-Up Histograph Parameters
numBins = 200;
maxEnergyBin = 700; % [keV]
binWidth = maxEnergyBin/numBins;
y1 = 0:binWidth:maxEnergyBin;

%[x0,y1] = hist(Etot,200);
[x1,y0] = hist(Edep,y1);
[x2,y2] = hist(Econtained, y1);
[x3,y3] = hist(Enoncontained, y1);
[x4,y4] = hist(E, y1);

% Save Histogram Info

sname = ['SummaryInfo', filenamePhrase];
save(sname, 'y0', 'y1', 'y2', 'y3', 'y4', 'x1', 'x2', 'x3', 'x4');

sname2 = ['SummaryEventInfo', filenamePhrase];
save(sname2, 'Etot', 'Edep', 'Eesc', 'Econtained','Enoncontained', 'multi');


% Plot Histogram Info
if plotFlag == 1
    %plot(y1, x1, 'g', y1,x2,'b', y1,x3,'r')
    
    % Calculate Initial Number of Gamma's
    nps = 101*500000 %4*10*10^6+2*5*10^6; %Total Number of Starting Photons
 
    
    figure(2)
    subplot(2,1,2);plot(y1, x1./nps, 'k', ...
         y1, x2./nps, 'b', ...
         y1, x3./nps, 'r','LineWidth',2);

    legend('Total Energy Spectrum','Cointained Electrons','Non-contained Electrons');
    title('Simulated CCD Energy Spectrum of Benchmarking Experiment', 'Fontsize',14)
    xlabel('CCD Event Energy [keV]','Fontsize',12);
    ylabel('Coincident Event Probability','Fontsize',12);
    %axes('Fontsize',20)
    nps2 = length(E);
    subplot(2,1,1);plot(y1,x4./nps2,'b','LineWidth',2);
    title('Experimental CCD Total Energy Spectrum','Fontsize',14);  
    xlabel('CCD Event Energy [keV]','Fontsize',12);
    ylabel('Event Probability [arb units]','Fontsize',12);

end