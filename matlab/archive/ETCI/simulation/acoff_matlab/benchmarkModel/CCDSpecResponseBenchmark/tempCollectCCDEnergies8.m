function res = tempCollectCCDEnergies8(Event)
%[Etot2, EimgTot2, Edep2, EescOld2, EescNew2, multi2, error1,error2] = collectCCDEnergies2(electronTRKArray)
% [Etot2, Edep2, EescOld2, EescNew2, multi2] =
% collectCCDEnergies2(electronTRKArray)
% errTest = calculateARM(electronTRKArray)
% Current as of: 3/14/12

%load(electronTRKArray);

% Pre-allocate Space for Saving CCD Events of Interest:
% Find Number of Energy trackable CCD Events
numCCDInteractions =[];
realTest =[];
res.error1 = [];

for i = 1:length(Event)
    try
    numCCDInteractions = [numCCDInteractions length(Event{1,i}.T)];
    realTest = [realTest isreal(Event{1,i}.E)];
    catch e
        res.error1 = [res.error1 e];
        continue 
    end
end

% Pre-Allocate Array Space for CCD Array
numCCD = sum(numCCDInteractions);  % length(Event); %
res = {};
    %res.Etot  = NaN(1,numCCD);
    res.EimgTot = NaN(1,numCCD);
    res.Edep  = NaN(1,numCCD);
    res.E0 = NaN(1,numCCD);
    res.Eesc = NaN(1,numCCD);
    res.Emulti = NaN(1,numCCD);
    res.multi = NaN(1,numCCD);
    res.EimgContained = NaN(1,numCCD);
    res.EimgNonContained = NaN(1,numCCD);
%     res.Ediff_G4_Diffusion = NaN(1,numCCD);
%     res.Ediff_G4_image = NaN(1,numCCD);
%     res.eventNumImg = NaN(1,numCCD);
%     res.eventNumDiffusion = NaN(1,numCCD);
%     res.EtotalGe = NaN(1,numCCD);
    n = 1;
    %
    
%     % Testing of DiffuseTrack Accurately reading G4 Tracking Matrix
%     % Etot - Edep: (G4 Electron Initial Energy) - (DiffuseTrack Edep) 
%     res.E_G4diff_Contained_Diffusion =NaN(1,numCCD);
%     res.E_G4diff_NonContained_Diffusion= NaN(1,numCCD);
%     % Etot - Eimg: (G4 Electron Initial Energy) - (CCDSegment Sum Image Energy)
%     res.E_G4diff_Contained_Image = NaN(1,numCCD);
%     res.E_G4diff_nonContatined_Image= NaN(1,numCCD);
    
    res.error1 = [];
    res.error2 = [];
    res.error3 = [];
    res.error4 = [];
    %numError4 = 0;
%     res.error4 = [];
%     res.error5 = [];
%     res.error6 = [];

for i = 1:numCCD %length(Event)
    
%     if isempty(Event{1,i}.out.cheat.Etot)   % Make sure the track energies are recorded
%         errTest.err(i) = 0;       % Zero Array if failed
%         contiune 
%     else
%         errTest.err(i) = 1;
%     end


% Collect All Information for Histograming: 
%
% Save Image Energy, if there is no energy, it gets skipped. 
    try
        if  ~isempty(Event{1,i}.E(1)) % Make sure image energy exits
            
            % Re-Add for accounting for Multiple Track Images (MS)
            if length(Event{1,i}.T) == 1
                % Image Total Energy 
                %res.EimgTot(1,n) = Event{1,i}.out.E(1);
                res.EimgTot(1,n) = Event{1, i}.T{1, 1}.E(1);
                % Tracked Electron - Perfect Electron Energy 
                res.Edep(1,n) = Event{1,i}.cheat.Edep(1);
                % Exact Electron Starting Energy
                res.E0(1,n) = Event{1,i}.cheat.E0(1);
                % One measure of Multiplicity 
                res.Emulti = Event{1, i}.T{1, 1}.E(1);
                res.multi = length(Event{1,i}.T);
                
                 % Sort Image Energies Into Contained and Non-Contained
                 %if exist('Event{1,i}.out.E','var')
                 res.Eesc(1,n) = (Event{1,i}.cheat.E0)-(Event{1,i}.cheat.Edep(1)); 
                 if (res.Eesc(1,n)) >  100 % [eV]  ... Don't want to count atomic photon escape 
                     res.EimgNonContained(1,n) = Event{1, i}.T{1, 1}.E(1);
                     %res.E_G4diff_nonContatined_Image(1,n)= (Event{1,i}.out.cheat.Etot - Event{1,i}.out.E(1));
                 else
                     res.EimgContained(1,n) = Event{1, i}.T{1, 1}.E(1);
                     %res.E_G4diff_Contained_Image(1,n) = (Event{1,i}.out.cheat.Etot - Event{1,i}.out.E(1));
                 end
                 n = n+1;
            else
                for j = 1:length(Event{1,i}.T) 
                    try
                        % Image Total Energy 
                        res.EimgTot(1,n) = Event{1, i}.T{1, j}.E(1);

                        % Sort Image Energies Into Contained and Non-Contained
                        res.Eesc(1,n) = (Event{1,i}.cheat.E0)-(Event{1,i}.cheat.Edep(1)); % Not the best flag for multiple scatter... 
                        res.Emulti = Event{1, i}.T{1, j}.E(1);
                        
                         if (res.Eesc(1,n)) >  100 % [eV]  ... Don't want to count atomic photon escape 
                             res.EimgNonContained(1,n) = Event{1, i}.T{1, j}.E(1);
                             %res.E_G4diff_nonContatined_Image(1,n)= (Event{1,i}.out.cheat.Etot - Event{1,i}.out.E(1));
                         else
                             res.EimgContained(1,n) = Event{1, i}.T{1, j}.E(1);
                             %res.E_G4diff_Contained_Image(1,n) = (Event{1,i}.out.cheat.Etot - Event{1,i}.out.E(1));
                         end
                        % Tracked Electron - Perfect Electron Energy 
                        n = n+1;
                    catch e
                        %res.error4 = [res.error4 length(Event{1,i}.cheat)];
                        res.error4 = [res.error4 e];
                    end
                        
                end
                
            end
        end
        
    catch e
        % Save Event that doesn't contain info. of Interest
        %sname = ['ErrorInfo_',num2str(i),'_' electronTRKArray];
        %ErrorEvent = Event{1,i};
        %save(sname, 'ErrorEvent');
         
        res.error2 = [res.error2 e];
        res.error3 = [res.error3 i];
%         %res.error3 = [res.error3 length(Event{1,i}.T)];
%         if  ~isempty(Event{1,i}.T)
%             %res.error4 = [res.error4 length(Event{1,i}.cheat)];
%             res.error4 = [res.error3 length(Event{1,i}.T)];
%         end
        
        %numError4 = numError4+1;
        %res.error4{1,numError4}.lengthT = length(Event{1,i}.T);
        %res.error4{1,numError4}.lengthCheat = length(Event{1,i}.cheat);
        %length(Event{1,i}.T)
        %length(Event{1,i}.cheat)
        % Temp Fix - Place a 0 placeholder... (?)
        %Don't want placeholder.. just skip the problem. 
%         EimgTot2 = [EimgTot2 0];
%         Etot2 = [Etot2 0];
%         Edep2 = [Edep2 0];
%         EescOld2 = [EescOld2 0];
%         EescNew2 = [EescNew2 0]; 
%         multi2 = [multi2 0];
    end

%   %%% Energy Comparisons Code %%% -- Could still be Usefull Testing Tool    
%     try
%         tempEdiff1 = (Event{1,i}.out.cheat.Etot - Event{1,i}.out.E);
%         res.Ediff_G4_image = [res.Ediff_G4_image tempEdiff1];
%         if tempEdiff1 < -5
%             res.eventNumImg = [res.eventNumImg i];
%         end  
%     catch
%         %
%     end
%     
%     try
%         res.Edep = [res.Edep Event{1,i}.out.cheat.Edep];
%         res.error3 = [res.error3 0];
%     catch 
%         res.error3 = [res.error3 1];
%     end
%     
%     try
%         tempEdiff2 = (Event{1,i}.out.cheat.Etot - Event{1,i}.out.cheat.Edep);
%         res.Ediff_G4_Diffusion = [res.Ediff_G4_Diffusion tempEdiff2];
%         if tempEdiff2 < -5
%             res.eventNumDiffusion = [res.eventNumDiffusion i];
%         end
%     catch
%         %
%     end
%     
%     try 
%         res.EescOld = [res.EescOld Event{1,i}.Eesc];  % Geant4 Track escape energy
%         res.error4 = [res.error4 0];
%     catch 
%         res.error4 = [res.error4 1];
%     end
%     
%     try
%         res.EescNew = [res.EescNew (Event{1,i}.out.cheat.Etot-Event{1,i}.out.cheat.Edep)]; % Not looking at this yet. 
%         res.error5 = [res.error5 0];
%     catch 
%         res.error5 = [res.error5 1];
%     end
%     
%     try
%        res.multi = [res.multi Event{1,i}.multiplicity];
%        res.error6 = [res.error6 0];
%     catch 
%         res.error6 = [res.error6 1];
%     end
    
end



