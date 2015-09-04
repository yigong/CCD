function out = collectEvents6(filenamePhrase,E, plotFlag, diagnosisFlag)
% [Etot, Edep, EescOld, EescNew, Econtained,Enoncontained,
% multi,error1,error2] = collectEvents2(filenamePhrase,E, plotFlag)
% [Etot, Edep, Eesc, Econtained,Enoncontained, multi] =
% collectEvents2(filenamePhrase,E) 
% filename == _TRK.mat file with Event{} Structure.%phrase characteristic
% of filename
% E = CCD Benchmarking Experiment Tracks Energy Matrix
% Coincident Data Set: E = Compiled Results of Coincident experiemnt
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
out.err = {};
primaryPhotons = 10^6; % Number of primary photons for each file

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
% ARM Info
DA = [];
E_CCD_Plus_Ge = [];
%E_Ge_DE = [];
EtotalGe = [];
E_ARM_CCD = [];


% Diagostics Values: 
Ediff_G4_image = [];
Ediff_G4_Diffusion = [];

% Testing of DiffuseTrack Accurately reading G4 Tracking Matrix
% Etot - Edep: (G4 Electron Initial Energy) - (DiffuseTrack Edep) 
E_G4diff_Contained_Diffusion =[];
E_G4diff_NonContained_Diffusion= [];
% Etot - Eimg: (G4 Electron Initial Energy) - (CCDSegment Sum Image Energy)
E_G4diff_Contained_Image = [];
E_G4diff_nonContatined_Image= [];
eventNumImg = [];
eventNumDiffusion = [];
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
        
        %ARM
        DA = [DA  res.DA];
        E_CCD_Plus_Ge = [E_CCD_Plus_Ge res.E_CCD_Plus_Ge];
        %E_GeDet_dE = [E_GeDet_dE res.E_GeDet_dE];
        EtotalGe = [EtotalGe res.EtotalGe];
        E_ARM_CCD = [E_ARM_CCD res.E_ARM_CCD];
        
        %Diagnostics
        Ediff_G4_image = [Ediff_G4_image res.Ediff_G4_image];
        Ediff_G4_Diffusion = [Ediff_G4_Diffusion res.Ediff_G4_Diffusion];
        
        E_G4diff_Contained_Diffusion = [E_G4diff_Contained_Diffusion res.E_G4diff_Contained_Diffusion];
        E_G4diff_NonContained_Diffusion = [E_G4diff_NonContained_Diffusion res.E_G4diff_NonContained_Diffusion];
        
        E_G4diff_Contained_Image = [E_G4diff_Contained_Image res.E_G4diff_Contained_Image];
        E_G4diff_nonContatined_Image = [E_G4diff_nonContatined_Image res.E_G4diff_nonContatined_Image];
        %Forgot...
        eventNumImg = [eventNumImg res.eventNumImg];
        eventNumDiffusion = [eventNumDiffusion res.eventNumDiffusion];
        
        %Recording of Try/catch Errors from collectCCDEnergies5.m
        error1 = [error1 res.error1];
        error2 = [error2 res.error2];
        error3 = [error3 res.error3];
        error4 = [error4 res.error4];
        error5 = [error5 res.error5];
        error6 = [error6 res.error6];
    
        disp(['Done Analyzing ',fnames{i},'...'])
        n = n+1;
        toc
        
    catch e
        disp(['Skipping Analysis of ', fnames{i}, '... Does not contain proper Event Structure... ']);
        toc
        out.err{i} = e;
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

% Look at difference in size of Events Recorded and Events Properly
% Diffused
length(Etot);
length(EimgTot);

% Energy Window the CCD and the Ge-Detector
windowed_E_ARM_CCD = [];
windowed_DA = [];
windowed_EsystemTotal = [];

try
   for i=1:(length(E_ARM_CCD))
       if ((E_ARM_CCD(i) >= 140) && (E_ARM_CCD(i) <= 360))
           if((EtotalGe(i)>=300) && (EtotalGe(i)<=520)) 
                windowed_E_ARM_CCD = [windowed_E_ARM_CCD E_ARM_CCD(i)];
                windowed_DA = [windowed_DA DA(i)];
                windowed_EsystemTotal = [windowed_EsystemTotal E_CCD_Plus_Ge(i)];
           end
       end
   end
   out.windowed_DA = windowed_DA;
   out.windowed_EsystemTotal = windowed_EsystemTotal;
   out.windowed_E_ARM_CCD = windowed_E_ARM_CCD;
catch e 
        out.err{i} = e;
end

 
% Set-Up Histograph Parameters
numBins = 300;
maxEnergyBin = 680; % [keV]
binWidth = maxEnergyBin/numBins;
y1 = 0:binWidth:maxEnergyBin;

%[x0,y1] = hist(Etot,200);
%[x1,y0] = hist(Edep,y1);
[x1,y0] = hist(windowed_E_ARM_CCD,y1);
[x2,y2] = hist(EimgContained, y1);
[x3,y3] = hist(EimgNonContained, y1);
%[x4,y4] = hist(E.E_CCD, y1);
%[xARM, yARM] = hist(windowed_DA);

% Save Histogram Info


sname = ['SummaryInfo', filenamePhrase];
save(sname, 'y0', 'y1', 'y2', 'y3', 'x1', 'x2', 'x3');

sname2 = ['SummaryEventInfo', filenamePhrase];
save(sname2, 'Etot', 'EimgTot','Edep', 'EescOld','EescNew', 'EimgContained','EimgNonContained','Ediff_G4_image', 'Ediff_G4_Diffusion', 'multi', 'nps', 'eventNumImg', 'eventNumDiffusion', 'DA', 'E_CCD_Plus_Ge', 'EtotalGe', 'E_ARM_CCD', 'windowed_DA','windowed_EsystemTotal');

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
%out.nps2 = length(E.E_CCD);
out.eventNumImg = eventNumImg;
out.eventNumDiffusion = eventNumDiffusion;
out.gARM = DA;
out.E_ARM_CCD = E_ARM_CCD;
out.EtotalGe = EtotalGe;
out.EsystemTotal = E_CCD_Plus_Ge;

% Plot Histogram Info
if plotFlag == 1
    %plot(y1, x1, 'g', y1,x2,'b', y1,x3,'r')
    
    % Calculate Initial Number of Gamma's
    %nps = 1; %100*10*10^6 %4*10*10^6+2*5*10^6; %Total Number of Starting Photons

    figure(3)
    subplot(3,1,1);plot(y1, x1./nps, 'k', ...
         y1, x2./nps, 'b', ...
         y1, x3./nps, 'r','LineWidth',2);

    legend('Total Windowed Energy Spectrum','Cointained Electrons','Non-contained Electrons', 'Location','NorthEast');
    title('Simulated CCD Energy Spectrum of Benchmarking Experiment', 'Fontsize',14)
    xlabel('CCD Event Energy [keV]','Fontsize',12);
    ylabel('Coincident Event Probability','Fontsize',12);
    %axes('Fontsize',20)
    nps2 = length(E.E_CCD);
    subplot(3,1,2);
    plot(y1,x4./nps2,'b', ...
         y1, x1./nps2, 'r', 'LineWidth',2);
         %y1, x1./length(Edep), 'r', 'LineWidth',2); %Normalize Ba-133 Spec
         % with total number of events in spectrum. 
    legend('Experimental Energy Spectrum', 'Simulated CCD Energy Spectrum', 'Location','NorthEast');
    title('Experimental CCD Total Energy Spectrum','Fontsize',14);  
    xlabel('CCD Event Energy [keV]','Fontsize',12);
    ylabel('Event Probability [arb units]','Fontsize',12);
    
    subplot(3,1,3);
    
    %plot(xARM,yARM, 'r')
    
    
    toc

   
    figure(6)
    %plot(E_CCD_Plus_Ge, DA, '.b');
    plot(E.Etot,E.dtheta2,'.r', windowed_EsystemTotal, windowed_DA , '.b'); 
    title('Simulated CCD+GeDSSD Coincident-Experiment Feature Comparison     ', 'Fontsize',14)
    xlabel('Total Energy Detected [keV]        ', 'Fontsize',12)
    ylabel('Angular Resolusion Measure (ARM) [deg]         ','Fontsize', 12);
    %axis([630 (662+32) -50 50 ]);
end
    

if diagnosisFlag
    tic
    figure(4)
    subplot(2,1,1);
    hist(Ediff_G4_Diffusion, -5:.01:5)
    title('Energy Difference between Geant4 Track(Etot) and Edep (Diffusion track) Etot-Edep')
    subplot(2,1,2)
    hist(Ediff_G4_image,-5:.01:5)
    title('Energy Difference between Geant4 Track(Etot) and Track Image(Eimg) Etot-Eimg')
    toc
    tic
    figure(5)
    % 4 Histograms in one plot... 
    subplot(2,2,1);
    hist(E_G4diff_Contained_Diffusion, -5:.01:5)
    title('E G4diff Contained Diffusion')
    
    subplot(2,2,3);
    hist(E_G4diff_NonContained_Diffusion,-400:.01:500)
    title('E G4diff NonContained Diffusion')
    
    subplot(2,2,2);
    hist(E_G4diff_Contained_Image,-5:.01:5)
    title('E G4diff Contained-Image')
    
    subplot(2,2,4)
    hist(E_G4diff_nonContatined_Image,-400:.01:500)
    title('E G4diff NonContatined - Image')
    toc
end
