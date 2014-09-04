function makeMacroFile2(macroName, folderLocation, saveName, startNum, endNum, sourceName, gunEnergyKeV, nps)
% makeMacroFile2(macroName, folderLocation, saveName, startNum, endNum,sourceName, gunEnergyKeV, nps)
% Examples and past inputs:
% makeMacroFile2('Ba133_ang1_1000nm_10nmRC-CFS_100mil_1-50_Batch_01', '/global/lr3/acoffer/g4Data/CCD_Tracking/Ba133_angle1', 'Ba133_ang1_1000nmSM_10nmRC-CFS_100mil', 1, 50,'Ba133', 0, 100000000)
% Batch02:  (10mil) /global/lr3/acoffer/g4Data/CCD_Tracking/Ba133_angle1/Ba133_ang1_10mil_Batch02
% makeMacroFile2('Ba133_ang1_10mil_1-50_Batch_01','/global/lr3/acoffer/g4Data/CCD_Tracking/Ba133_angle1/Ba133_ang1_10mil_Batch02','Ba133_ang1_10nmRC-CFS_10mil_1-50_Batch_01',1,50,'Ba133',0,10000000);
% New build Ba-133_ang1 (100nm RC-CFS) (10mil batches) 
% /global/lr3/acoffer/g4Data/CCD_Tracking/Ba133_angle1/Ba133_ang1_100nmRC-CFS_10mil_Batch03
% makeMacroFile2('Ba133_ang1_100nmRC-CFS_10mil_1-50_Batch_01','/global/lr3/acoffer/g4Data/CCD_Tracking/Ba133_angle1/Ba133_ang1_100nmRC-CFS_10mil_Batch03','Ba133_ang1_100nmRC-CFS_10mil_1-50',1,50,'Ba133',0,10000000);
fileSaveName = [macroName '.mac'];
if isempty(who('gunEnergyKeV'))
    gunEnergyKeV = 0; % ion sources need a gps gun energy of 0
end
if strcmpi(sourceName, 'Ba133')
    sourceVolumePosition = '/testem/det/setSourcePosition -187.07 29.02 140.85 mm \n';
    sourcePosition = '/gps/position -187.07 29.02 140.85 mm \n\n/grdm/decayDirection 1 0 0  #* Supply the direction vector for decay products\n/grdm/decayHalfAngle 180 deg \n';
    radDecay = true;
else
    sourceVolumePosition = '/testem/det/setSourcePosition -20.0 0.0 0.0 mm \n';
    sourcePosition = '/gps/position -925.2782 -8.3933 374.5824 mm \n/gps/direction 0.92718 0.0041723 -0.37458 \n'; % Pencil Beam start pt. 
    radDecay = false;
end

%
% Non-Repeated Information:
fid = fopen(fileSaveName,'w');
    % Header Info
    fprintf(fid,'# ');
    fprintf(fid, fileSaveName);
    fprintf(fid,'\n\n/control/verbose 2 \n/run/verbose 2 \n');
    fprintf(fid,sourceVolumePosition);
    fprintf(fid,'\n/run/initialize \n');
   % StepMax Settings:
    fprintf(fid, '\n# Set StepMax in Calrimeter Layers:\n');
    fprintf(fid,'/testem/stepMax/absorber 1 1000 nm \n');
    fprintf(fid,'/testem/stepMax/absorber 2 1000 nm \n');
    fprintf(fid,'/testem/stepMax/absorber 3 1000 nm \n');
    
    %Source Volume
    % GPS Source Postion
    % Run initialize 
    if radDecay
        % Source
        fprintf(fid,'\n/gps/particle ion \n# Format:  /gps/ion <Z A Q E> \n/gps/energy 0. MeV \n/gps/ion 56 133 0 0. \n');
        fprintf(fid, sourcePosition);
    else
        fprintf(fid,'/gps/particle gamma \n/gps/energy 661.657 keV \n');
    end

    % Filename Location:
    fprintf(fid,'\n/SimulationOutput/FileLocation ');
    fullFolderLocation = [folderLocation '/'];
    fprintf(fid, fullFolderLocation);
    fprintf(fid, '\n/testem/event/printModulo 500000 \n');
    
% Run Information for each nps per run
for  i=(startNum:endNum)
    fprintf(fid,'\n/SimulationOutput/FileName ');
    binaryFileName = [saveName '_0' num2str(i) '.dat'];
    fprintf(fid, binaryFileName);
    fprintf(fid, '\n/run/beamOn ');
    numParts = num2str(nps);
    fprintf(fid, numParts);
    fprintf(fid, '\n');

end

status = fclose(fid);
clear fid;

%  %New Write to File Commands look like: BA133
% 
% /control/verbose 2
% /run/verbose 2
% #
% /testem/det/setSourcePosition -187.07 29.02 140.85 mm
% /run/initialize
% 
% # Set StepMax in Calrimeter Layers
% /testem/stepMax/absorber 1 1000 nm
% /testem/stepMax/absorber 2 1000 nm
% /testem/stepMax/absorber 3 1000 nm
% 
% # SOURCE
% /gps/particle ion
% # Format:  /gps/ion <Z A Q E>
% /gps/energy 0. MeV
% /gps/ion 56 133 0 0.
% /gps/position -187.07 29.02 140.85 mm
% /gps/direction 1 0 0
% 
% /grdm/decayDirection 1 0 0  #* Supply the direction vector for decay products
% /grdm/decayHalfAngle 180 deg
% 
% /SimulationOutput/FileLocation /global/lr3/acoffer/g4Data/CCD_Tracking/Ba133_angle1/
% /SimulationOutput/FileName gpsBa133_10mil_01.dat
% #
% /testem/event/printModulo 1000000
% #
% #/run/beamOn 100
% # 100 Million
% /run/beamOn 100000000