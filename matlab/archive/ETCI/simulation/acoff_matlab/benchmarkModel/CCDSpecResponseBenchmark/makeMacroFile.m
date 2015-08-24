function makeMacroFile(macroName, folderLocation, saveName, startNum, endNum, nps)
% makeMacroFile(macroName, folderLocation, saveName, startNum, endNum, nps)
fileSaveName = [macroName '.mac'];

fid = fopen(fileSaveName,'w');
    % Header Info
    fprintf(fid,'/control/verbose 2 \n/control/saveHistory \n/run/verbose 2 \n/tracking/verbose 0\n');
    

% Run Information for each 100000 run
for  i=(startNum:endNum)
    fprintf(fid, '\n/gun/energy .662 MeV');
    fprintf(fid, '\n/run/initialize \n/process/eLoss/fluct 0 \n/SimulationOutput/FileName ' );
    fprintf(fid, folderLocation);
    binaryFileName = [saveName num2str(i) '.dat'];
    fprintf(fid, binaryFileName);
    fprintf(fid, '\n/run/beamOn ');
    numParts = num2str(nps);
    fprintf(fid, numParts);
    fprintf(fid, '\n');

end

status = fclose(fid);
clear fid;


% /gun/energy .662 MeV
% /step/file/name /home/data2/eTCCI_Data_Feb_2010/binaryBackgroundData/eT_Gam_dirt_1K_test3.bin
% /step/file/clear
% /run/initialize
% /process/eLoss/fluct 0
% /run/beamOn 1000
% /step/file/write
% 
% /gun/energy .662 MeV
% /step/file/name /home/data2/eTCCI_Data_Feb_2010/binaryBackgroundData/eT_Gam_dirt_10K_test4.bin
% /step/file/clear
% /run/initialize
% /process/eLoss/fluct 0
% /run/beamOn 10000
% /step/file/write
% 
% /gun/energy .662 MeV
% /step/file/name /home/data2/eTCCI_Data_Feb_2010/binaryBackgroundData/eT_Gam_dirt_1.bin
% /step/file/clear
% /run/initialize
% /process/eLoss/fluct 0
% /run/beamOn 100000
% /step/file/write
% 
% /gun/energy .662 MeV
% /step/file/name /home/data2/eTCCI_Data_Feb_2010/binaryBackgroundData/eT_Gam_dirt_10.bin
% /step/file/clear
% /run/initialize
% /process/eLoss/fluct 0
% /run/beamOn 100000
% /step/file/write
% 
% /gun/energy .662 MeV
% /step/file/name /home/data2/eTCCI_Data_Feb_2010/binaryBackgroundData/eT_Gam_dirt_11.bin
% /step/file/clear
% /run/initialize
% /process/eLoss/fluct 0
% /run/beamOn 100000
% /step/file/write
% 
% /gun/energy .662 MeV
% /step/file/name /home/data2/eTCCI_Data_Feb_2010/binaryBackgroundData/eT_Gam_dirt_12.bin
% /step/file/clear
% /run/initialize
% /process/eLoss/fluct 0
% /run/beamOn 100000
% /step/file/write


