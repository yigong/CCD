function out = collectEvents5(filenamePhrase,E, plotFlag)
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
% Current as of: 3/14/12
                                            
IDphrase = filenamePhrase;
f = dir;
n=1;
%Event = {};
out = {};
primaryPhotons = 10^7; % Number of primary photons for each file

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
Ediff_G4_image = [];
Ediff_G4_Diffusion = [];
% Testing of DiffuseTrack Accurately reading G4 Tracking Matrix
% Etot - Edep: (G4 Electron Initial Energy) - (DiffuseTrack Edep) 
E_G4diff_Contained_Diffusion =[];
E_G4diff_NonContained_Diffusion= [];
% Etot - Eimg: (G4 Electron Initial Energy) - (CCDSegment Sum Image Energy)
E_G4diff_Contained_Image = [];
E_G4diff_nonContatined_Image= [];

    error1 = [];
    error2 = [];
    error3 = [];
    error4 = [];
    error5 = [];
    error6 = [];
    
n= 0;

for i=1:length(fnames)
    disp(['Analyzing',fnames{i},'...'])
    tic
    try
        
        %[tempEtot,tempEimgTot, tempEdep,tempEescOld, tempEescNew, tempMulti, tempError1,tempError2] = collectCCDEnergies2(fnames{i}); % Saves calculated events
        res = collectCCDEnergies5(fnames{i}); % Saves calculated events
     
        Etot = [Etot res.Etot];
        EimgTot = [EimgTot res.EimgTot];
        Edep = [Edep res.Edep];
        
        EescOld = [EescOld res.EescOld];
        EescNew= [EescNew res.EescNew];
        multi = [multi res.multi];
        
        EimgNonContained = [EimgNonContained res.EimgNonContained];
        EimgContained = [EimgContained res.EimgContained];
        
        Ediff_G4_image = [Ediff_G4_image res.Ediff_G4_image];
        Ediff_G4_Diffusion = [Ediff_G4_Diffusion res.Ediff_G4_Diffusion];
        
        E_G4diff_Contained_Diffusion = [E_G4diff_Contained_Diffusion res.E_G4diff_Contained_Diffusion];
        E_G4diff_NonContained_Diffusion = [E_G4diff_NonContained_Diffusion res.E_G4diff_NonContained_Diffusion];
        
        E_G4diff_Contained_Image = [E_G4diff_Contained_Image res.E_G4diff_Contained_Image];
        E_G4diff_nonContatined_Image = [E_G4diff_nonContatined_Image res.E_G4diff_nonContatined_Image];
        
        
        
        error1 = [error1 res.error1];
        error2 = [error2 res.error2];
        error3 = [error3 res.error3];
        error4 = [error4 res.error4];
        error5 = [error5 res.error5];
        error6 = [error6 res.error6];
       
        %tempEtot= []; tempEimgTot = []; tempEdep=[]; tempEescOld=[]; tempEescNew; tempMulti=[]; tempError1 = []; tempError2 = [];
    
        disp(['Done Analyzing ',fnames{i},'...'])
        n = n+1;
        toc
        
    catch
        disp(['Skipping Analysis of ', fnames{i}, '... Does not contain proper Event Structure... ']);
        toc
        continue
        
    end
end

% Calculate Number of Initial Photons
nps = n*primaryPhotons;
sumError1 = sum(error1)
sumError2 = sum(error2)
sumError3 = sum(error3)
sumError4 = sum(error4)
sumError5 = sum(error5)
sumerror6 = sum(error6)

if plotFlag ==1
    % Plot Energy Spectrum of Interest
    disp(['Plotting ', filenamePhrase, '...'])
    tic
end

% Determine Energy Spec. of Cointained and Non-Cointained Electrons
% Econtained = [];
% Enoncontained = [];
length(Etot)
length(EimgTot)

% Work in Progress.. 3/14/12... 
% Need to move Contained/ non-containe tracking to collectCCDEnergies
% Section...
%  for i=1:length(EimgTot)%length(Etot)
%      if abs(EescOld(i)) > 0     % Escape Determined by Geant4 Track and Etot vs Edep. 
%          % Enoncontained = [Enoncontained Edep(i)]; % Geant4 Track Energy
%          Enoncontained = [Enoncontained EimgTot(i)]; % DiffuseTrack Image Energy
%          
%      else
%          %Econtained = [Econtained Edep(i)];    % Geant4 Track Energy
%          Econtained = [Econtained EimgTot(i)];  % DiffuseTrack Image Energy
%      end
%  end
 
% Set-Up Histograph Parameters
numBins = 300;
maxEnergyBin = 680; % [keV]
binWidth = maxEnergyBin/numBins;
y1 = 0:binWidth:maxEnergyBin;

%[x0,y1] = hist(Etot,200);
%[x1,y0] = hist(Edep,y1);
[x1,y0] = hist(EimgTot,y1);
[x2,y2] = hist(EimgContained, y1);
[x3,y3] = hist(EimgNonContained, y1);
[x4,y4] = hist(E, y1);

% Save Histogram Info


sname = ['SummaryInfo', filenamePhrase];
save(sname, 'y0', 'y1', 'y2', 'y3', 'y4', 'x1', 'x2', 'x3', 'x4');

sname2 = ['SummaryEventInfo', filenamePhrase];
save(sname2, 'Etot', 'EimgTot','Edep', 'EescOld','EescNew', 'EimgContained','EimgNonContained','Ediff_G4_image', 'Ediff_G4_Diffusion', 'multi', 'nps');

out.Etot = Etot;
out.EimgTot = EimgTot;
out.Edep = Edep;
out.EescOld = EescOld;
out.EescNew = EescNew;
out.EimgContained = EimgContained;
out.EimgNonContained = EimgNonContained;
out.Ediff_G4_Diffusion = Ediff_G4_Diffusion;
out.Ediff_G4_image = Ediff_G4_image;
out.multi = multi;
out.errorType1 = error1;
out.errorType2 = error2;
out.nps = nps;

% Plot Histogram Info
if plotFlag == 1
    %plot(y1, x1, 'g', y1,x2,'b', y1,x3,'r')
    
    % Calculate Initial Number of Gamma's
    %nps = 1; %100*10*10^6 %4*10*10^6+2*5*10^6; %Total Number of Starting Photons
 
    
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
    tic
    figure(4)
    subplot(2,1,1);
    hist(Ediff_G4_Diffusion, -5:.01:5)
    title('Energy Difference between Geant4 Track(Etot) and Edep (Diffusion track) Etot-Edep')
    subplot(2,1,2)
    hist(Ediff_G4_image,-5:.01:5)
    title('Energy Difference between Geant4 Track(Etot) and Track Image(Eimg) Etot-Eimg')
    toc
    
end

