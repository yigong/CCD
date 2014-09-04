function res = collectCCDEnergies5(electronTRKArray)
%[Etot2, EimgTot2, Edep2, EescOld2, EescNew2, multi2, error1,error2] = collectCCDEnergies2(electronTRKArray)
% [Etot2, Edep2, EescOld2, EescNew2, multi2] =
% collectCCDEnergies2(electronTRKArray)
% errTest = calculateARM(electronTRKArray)
% Current as of: 3/14/12

load(electronTRKArray);
n = length(Event);

res = {};
    res.Etot  = [];
    res.EimgTot = [];
    res.Edep  = [];
    res.EescOld  = [];
    res.EescNew = [];
    res.multi = [];
    res.EimgContained = [];
    res.EimgNonContained = [];
    res.Ediff_G4_Diffusion = [];
    res.Ediff_G4_image = [];
    res.eventNumImg = [];
    res.eventNumDiffusion = [];
    res.EtotalGe = [];
    
    %
    % Purly ARM Output Set
    % ARMinfo.theta = [];
    % ARMinfo.gARM = [];
    % ARMinfo.gARMdeg = [];
    % ARMinfo.EtotalSystem = [];
    % ARMinfo.EtotalGe = [];
    % ARMinfo.E_ARM_CCD = [];

    res.DA = []; % Photon ARM Info (ARMinfo.gARMdeg)
    res.E_CCD_Plus_Ge = []; % For dA/Etot plots (ARMinfo.EtotalSystem )
    res.E_ARM_CCD = []; % CCD Energy Used for ARM (ARMinfo.E_ARM_CCD)
    
    % Testing of DiffuseTrack Accurately reading G4 Tracking Matrix
    % Etot - Edep: (G4 Electron Initial Energy) - (DiffuseTrack Edep) 
    res.E_G4diff_Contained_Diffusion =[];
    res.E_G4diff_NonContained_Diffusion= [];
    % Etot - Eimg: (G4 Electron Initial Energy) - (CCDSegment Sum Image Energy)
    res.E_G4diff_Contained_Image = [];
    res.E_G4diff_nonContatined_Image= [];
    
    res.error1 = [];
    res.error2 = [];
    res.error3 = [];
    res.error4 = [];
    res.error5 = [];
    res.error6 = [];

for i = 1:n
    
%     if isempty(Event{1,i}.out.cheat.Etot)   % Make sure the track energies are recorded
%         errTest.err(i) = 0;       % Zero Array if failed
%         contiune 
%     else
%         errTest.err(i) = 1;
%     end

% Calculate things of Interest that Haven't Been Saved in Event{}
%


% Collect All Information for Histograming: 
%
% Save Image Energy, if there; is no energy, it gets skipped. 
    try
         res.EimgTot = [res.EimgTot Event{1,i}.out.E];
         
         % Calculate ARM Info:
         % ARM Calculation only works if DiffuseTrack5.m returns an Electron
         % CCD Track Image energy: Event{1,i}.out.E 
         ARMinfo = calculateARM5(Event{1,i});
         
         res.DA = [res.DA ARMinfo.gARMdeg]; % Photon ARM Info
         res.E_CCD_Plus_Ge = [res.E_CCD_Plus_Ge ARMinfo.E_CCD_Plus_Ge]; % For dA/Etot plots
         %res.E_GeDet_dE = [res.E_GeDet_dE ARMinfo.E_GeDet_dE];
         res.EtotalGe =[res.EtotalGe ARMinfo.EtotalGe];
         res.E_ARM_CCD = [res.E_ARM_CCD ARMinfo.E_ARM_CCD];
         
         % Sort Image Energies Into Contained and Non-Contained
         %if exist('Event{1,i}.out.E','var')
         if (Event{1,i}.Eesc) > 0
             res.EimgNonContained = [res.EimgNonContained Event{1,i}.out.E];
             res.E_G4diff_nonContatined_Image= [res.Ediff_G4_image (Event{1,i}.out.cheat.Etot - Event{1,i}.out.E)];
         else
             res.EimgContained = [res.EimgContained Event{1,i}.out.E];
             res.E_G4diff_Contained_Image = [res.Ediff_G4_image (Event{1,i}.out.cheat.Etot - Event{1,i}.out.E)];
         end
         %end
         res.error1 = [res.error1 0];
    catch 
        res.error1 = [res.error1 1];
    end
    
    try
         % Sort Diffuse Edep Energies Into Contained and Non-Contained

         if (Event{1,i}.Eesc) > 0
             res.E_G4diff_NonContained_Diffusion =  [res.Ediff_G4_Diffusion (Event{1,i}.out.cheat.Etot - Event{1,i}.out.cheat.Edep)];
         else
             res.E_G4diff_Contained_Diffusion = [res.Ediff_G4_Diffusion (Event{1,i}.out.cheat.Etot - Event{1,i}.out.cheat.Edep)];
         end
         
    catch
        %
    end
    

    try
        res.Etot = [res.Etot Event{1,i}.out.cheat.Etot];
        % res.Edep = [res.Edep Event{1,i}.out.cheat.Edep];
        %res.EescOld = [res.EescOld Event{1,i}.Eesc];  % Geant4 Track escape energy
        %res.EescNew = [res.EescNew (Event{1,i}.out.cheat.Etot-Event{1,i}.out.cheat.Edep)]; % Not looking at this yet. 
        %res.multi = [res.multi2 Event{1,i}.multiplicity];
        res.error = [res.error2 0];
    catch 
        % Save Event that doesn't contain info. of Interest
        %sname = ['ErrorInfo_',num2str(i),'_' electronTRKArray];
        %ErrorEvent = Event{1,i};
        %save(sname, 'ErrorEvent');
        res.error2 = [res.error2 1];
        
        % Temp Fix - Place a 0 placeholder... (?)
        %Don't want placeholder.. just skip the problem. 
%         EimgTot2 = [EimgTot2 0];
%         Etot2 = [Etot2 0];
%         Edep2 = [Edep2 0];
%         EescOld2 = [EescOld2 0];
%         EescNew2 = [EescNew2 0]; 
%         multi2 = [multi2 0];
    end
    
    try
        tempEdiff1 = (Event{1,i}.out.cheat.Etot - Event{1,i}.out.E);
        res.Ediff_G4_image = [res.Ediff_G4_image tempEdiff1];
        if tempEdiff1 < -5
            res.eventNumImg = [res.eventNumImg i];
        end  
    catch
        %
    end
    
    try
        res.Edep = [res.Edep Event{1,i}.out.cheat.Edep];
        res.error3 = [res.error3 0];
    catch 
        res.error3 = [res.error3 1];
    end
    
    try
        tempEdiff2 = (Event{1,i}.out.cheat.Etot - Event{1,i}.out.cheat.Edep);
        res.Ediff_G4_Diffusion = [res.Ediff_G4_Diffusion tempEdiff2];
        if tempEdiff2 < -5
            res.eventNumDiffusion = [res.eventNumDiffusion i];
        end
    catch
        %
    end
    
    try 
        res.EescOld = [res.EescOld Event{1,i}.Eesc];  % Geant4 Track escape energy
        res.error4 = [res.error4 0];
    catch 
        res.error4 = [res.error4 1];
    end
    
    try
        res.EescNew = [res.EescNew (Event{1,i}.out.cheat.Etot-Event{1,i}.out.cheat.Edep)]; % Not looking at this yet. 
        res.error5 = [res.error5 0];
    catch 
        res.error5 = [res.error5 1];
    end
    
    try
       res.multi = [res.multi Event{1,i}.multiplicity];
       res.error6 = [res.error6 0];
    catch 
        res.error6 = [res.error6 1];
    end
    
end



