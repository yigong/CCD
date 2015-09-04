function G4THworker(loadpath,flist,savepath,savepre,loadpre,loadsuf,savesuf,noise,pixelSize,pixelNoiseMatrix,i,placeholdersuf,psft)
%function G4THworker(loadpath,flist,savepath,savepre,loadpre,loadsuf,savesuf,noise,pixelSize,pixelNoiseMatrix,i,placeholdersuf,psft)
%
% called by G4THparallel.m
%
% revised Jan 3, 2013

fileload = load(fullfile(loadpath,flist(i).name));

addpath /global/scratch/bcplimley/multi_angle/MatlabFunctions/

%batches of 1000, then save off
batchsize = 1000;
    for j=1:ceil(length(fileload.Event)/batchsize)
        numstring = [flist(i).name(length(loadpre)+1:strfind(flist(i).name,loadsuf)-1),'_',num2str(j)];  %includes flist number as well as part number
        savename = [savepre,numstring,savesuf];
        placeholdername = [savepre,numstring,placeholdersuf];
        if ~isempty(dir(fullfile(savepath,savename)))
            disp(['Skipping ',savename,' at ',datestr(now)])
            continue
        elseif ~isempty(dir(fullfile(savepath,placeholdername)))
            disp(['Skipping ',placeholdername,' at ',datestr(now)])
            continue
        end
        
        %write placeholder file
        save(fullfile(savepath,placeholdername),'numstring');
        
        disp(['Starting ',savename,' at ',datestr(now),' with ',num2str(CheckMemUse(whos)),' bytes in memory'])
        
        %each event
        for k = (j-1)*batchsize+1:min(j*batchsize,length(fileload.Event))
%             filesave.CCD{k-(j-1)*batchsize} = fileload.Event{k};  %not anymore.
            try
                DT{k-(j-1)*batchsize} = ...
                    Geant4TrackHandling7(fileload.Event{k}.trackM,'coinc',...
                    'driftdim','z',...
                    'depthcoordinates',[-0.65,0],...
                    'noise',noise.*1e-3,...
                    'pixelsize',pixelSize,...
                    'pixelnoisematrix',pixelNoiseMatrix,...
                    'psft',psft);
                %don't keep all the redundant copies of cheat.
                fieldNames = fieldnames(DT{k-(j-1)*batchsize});
                %keep one full copy on top of the structure.
                currentCheat = DT{k-(j-1)*batchsize}.(fieldNames{1}).cheat;
                for m = 1:length(fieldNames)
                    if ~strcmp(fieldNames{m}(1:3),'pix')
                        %skip over fields that are not a pixelsize/noise entry.
                        continue
                    end
                    DT{k-(j-1)*batchsize}.(fieldNames{m}).cheat = ...
                        rmfield(DT{k-(j-1)*batchsize}.(fieldNames{m}).cheat, ...
                        {'x','dE','particleID','longStepLength','alphaLong','betaLong',...
                        'sourcePhotonE1','sourcePhotonE2','sourcePhotonDirection1','sourcePhotonDirection2'});
                end
                DT{k-(j-1)*batchsize}.cheat = currentCheat;
                %save other primary info here too. it is saved only once for all pixelsize/noise combinations.
                DT{k-(j-1)*batchsize}.trackM = fileload.Event{k}.trackM;
                DT{k-(j-1)*batchsize}.out = fileload.Event{k}.out;
                DT{k-(j-1)*batchsize}.Etot = fileload.Event{k}.Etot;
                DT{k-(j-1)*batchsize}.Edep = fileload.Event{k}.Edep;
                DT{k-(j-1)*batchsize}.Eesc = fileload.Event{k}.Eesc;
                DT{k-(j-1)*batchsize}.multiplicity = fileload.Event{k}.multiplicity;
                
                
            catch err
                DT{k-(j-1)*batchsize}.err = err;
                disp(['Warning: encountered error "',err.message,'" in ',flist(i).name,' event # ',num2str(k)])
            end
        end

        %save
%        disp(['Saving ',savename,' at ',datestr(now)])
        %previously used function "savefile" that i wrote... when G4THworker was still part of G4THparallel
%         savefile(fullfile(savepath,savename),filesave.DT,flist);
        save(fullfile(savepath,savename),'DT','flist')
        
        disp(['Finished ',savename,' at ',datestr(now),' with ',num2str(CheckMemUse(whos)),' bytes in memory'])
        %remove placeholder
        delete(fullfile(savepath,placeholdername));
    end

