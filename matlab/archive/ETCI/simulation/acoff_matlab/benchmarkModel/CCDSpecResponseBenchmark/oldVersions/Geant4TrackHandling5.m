function out = Geant4TrackHandling5(varargin)
%function out = Geant4TrackHandling(datFilename,mode)
%function out = Geant4TrackHandling(matFilename,mode)
%function out = Geant4TrackHandling(g4matrix,mode)
%function out = Geant4TrackHandling(...,argID,arg)
%
% Modified 4/13/2012. Fixed bug that skipped certain binding energy Edeps from multiple scattering of a single photon.
%                     Also cleaned up some oldcodeflag sections that could never be reached.
%                     Now known as Geant4TrackHandling5.m for compliance with Amy's numbering system.
% Modified 3/7/2012. May still have trouble loading dat or mat files with multiple events..
%
% Depends on:
%   DiffuseTrack5.m
%   AddCCDImageNoise.m
%   psft.mat, unless loaded as an optional argument
%
% datFilename: path/filename of a *.dat ASCII data file from geant4. (see format below)
% matFilename: path/filename of a *.mat Matlab data file with geant4 matrix in it. (see format below)
% g4matrix: numerical array of geant4 data. Should be 15 columns as follows:
%       Columns: (units)
%       1: Track ID
%       2: Parent ID
%       3: Step num
%       4: Charge                           (e-)
%       5,6,7: Initial position (x,y,z)     (mm)
%       8,9,10: Final position (x,y,z)      (mm)
%       11: Accumulated track length        (mm)
%       12: Step length                     (mm)
%       13: Final energy (after step)       (eV)
%       14: Initial energy (before step)    (eV)
%       15: Energy deposited (dE)           (eV)
%     default coordinate system is CCD depth along x dimension.
%
% mode: string to specify what type of simulation the g4matrix comes from.
%   'coinc': coincidence experiment
%           e.g.: coincBench_5mil_008.dat
%   'single': single CCD with photons incident upon it. (algorithm tests)
%           e.g.: eTracks_Si_1501.mat
% argID: string to identify optional argument
% arg: value of optional argument
%
% argID list:
%  psft: cell array of PSFs, or path/filename of psft.mat (defaults to psft.mat in current directory)
%  fullimage: set flag for saving the entire CCD image into the output structure 
%     (don't use in a loop, for memory's sake)
%  nonoise: set flag for saving CCD images without statistical or dark current noise
%  zcentered: set flag for using 0 as the center of the depth of the detector, rather than one face.
%  progressbar: set flag for using 'progressbar' function for multi-track operation
%  CCDsize: [x,y] size of the CCD in pixels. Default 726x1454 for coinc mode.
%     'var': set flag for adapting simulated size of device based on the event. 
%     default 'var' for 'single' mode.)
%  manualenergycorrection: set flag for manually correcting deposited energy 
%     from secondary electrons, to avoid double counting. (depends on geant4 step writing)
%  diffusetrack2b: set flag for using old diffusion code. (version 2b, which handles the secondary corrections)
%  noiselist: supply a vector of noise sigmas to use for output.
%  energywindow: supply [Emin,Emax] for events to process. units: keV (converted to eV for geant4 comarison)
%     The energy window applies to the initial energy of the primary electron of the source photon. Escaped energy is not subtracted.
%  depthdim: provide either {'x','y','z'} or {1,2,3} to define which dimension in geant4 is the thickness of the CCD.
%     default for coinc mode: z
%     default for single mode: x
%  negdepth: set flag for using [-0.650, 0] of depthdim for the CCD, rather than [0, +0.650].
%     default for coinc mode: true
%     default for single mode: false


%% Input handling

if nargin<2
    error('Need input data and mode as input arguments')
end

%data
if ischar(varargin{1}) && ~isempty(strfind(varargin{1},'.dat'))
    M = load(varargin{1},'-ascii');
elseif ischar(varargin{1}) && ~isempty(strfind(varargin{1},'.mat'))
    Mtemp = load(varargin{1});
    if length(fieldnames(Mtemp))==1
        tmp2 = fieldnames(Mtemp);
        tmp = Mtemp.(tmp2{1});
        M = tmp{1};
        clear Mtemp tmp
    else
        error('don''t know how to handle multiple variables in data file')
    end
elseif isnumeric(varargin{1})
    M = varargin{1};
end

%mode
if ~ischar(varargin{2})
    error('Mode should be a string')
elseif strcmpi(varargin{2},'coinc') || strcmpi(varargin{2},'coincidence')
    inputmode = 'coinc';
elseif strcmpi(varargin{2},'single') || strcmpi(varargin{2},'singleccd') || strcmpi(varargin{2},'single ccd')
    inputmode = 'single';
else
    error('Mode not recognized')
end

out = [];
psft = [];
fullimageflag = false;  %default
nonoiseflag = false;    %default is to add normal noise to image.
progressbarflag = false;
energyaccountingflag = false;
oldcodeflag = false;
zcenteredflag = false;
noiselist = [];
energywindow = [-Inf,Inf];

if strcmp(inputmode,'coinc')
    variablesizeflag = false;
elseif strcmp(inputmode,'single')
    variablesizeflag = true;
end
%optional arguments
i=3;
while true
    if i>nargin
        break
    elseif strcmpi(varargin{i},'psft')
        if nargin>i && ischar(varargin{i+1})
            tmp = load(varargin{i+1});
            psft = tmp.psft;
        elseif nargin>i && iscell(varargin{i+1})
            psft = varargin{i+1};
        elseif nargin==i
            warning('Psft input not found; looking in current directory')
        else
            error('Psft input must be cell array or filename')
        end
        i=i+2;
    elseif strcmpi(varargin{i},'fullimage') || strcmpi(varargin{i},'fullimg') || ...
            strcmpi(varargin{i},'fullimageflag') || strcmpi(varargin{i},'fullimgflag')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==1 || varargin{i+1}==0)))
            fullimageflag = logical(varargin{i+1});
            i=i+2;
        else
            fullimageflag = true;
            i=i+1;
        end
    elseif strcmpi(varargin{i},'nonoise') || strcmpi(varargin{i},'nonoiseflag')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==1 || varargin{i+1}==0)))
            nonoiseflag = logical(varargin{i+1});
            i=i+2;
        else
            nonoiseflag = true;
            i=i+1;
        end
    elseif strcmpi(varargin{i},'progressbar') || strcmpi(varargin{i},'progressbarflag')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==1 || varargin{i+1}==0)))
            progressbarflag = logical(varargin{i+1});
            i=i+2;
        else
            progressbarflag = true;
            i=i+1;
        end
    elseif strcmpi(varargin{i},'variablesize') || strcmpi(varargin{i},'variableCCDsize')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==1 || varargin{i+1}==0)))
            variablesizeflag = logical(varargin{i+1});
            i=i+2;
        else
            variablesizeflag = true;
            i=i+1;
        end
    elseif strcmpi(varargin{i},'CCDsize')
        if nargin>i && isnumeric(varargin{i+1}) && length(varargin{i+1})==2
            CCDsize = varargin{i+1};
            i=i+2;
        elseif nargin>i && ischar(varargin{i+1}) && (strcmpi(varargin{i+1},'var') || strcmpi(varargin{i+1},'variable'))
            variablesizeflag = true;
        end
    elseif strcmpi(varargin{i},'energyaccounting') || strcmpi(varargin{i},'manualenergycorrect') || ...
            strcmpi(varargin{i},'manualenergycorrectionflag') || strcmpi(varargin{i},'manualenergycorrection') || ...
            strcmpi(varargin{i},'energyaccountingflag')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==1 || varargin{i+1}==0)))
            energyaccountingflag = logical(varargin{i+1});
            i=i+2;
        else
            energyaccountingflag = true;
            i=i+1;
        end
    elseif strcmpi(varargin{i},'diffusetrack2b') || strcmpi(varargin{i},'oldcode')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==1 || varargin{i+1}==0)))
            oldcodeflag = logical(varargin{i+1});
            i=i+2;
        else
            oldcodeflag = true;
            i=i+1;
        end
    elseif strcmpi(varargin{i},'zcentered') || strcmpi(varargin{i},'zcenteredflag')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==1 || varargin{i+1}==0)))
            zcenteredflag = logical(varargin{i+1});
            i=i+2;
        else
            zcenteredflag = true;
            i=i+1;
        end
    elseif strcmpi(varargin{i},'noiselist')
        if nargin>i && isnumeric(varargin{i+1})
            noiselist = varargin{i+1};
            i=i+2;
        else
            error('problem with noiselist input argument')
        end
    elseif strcmpi(varargin{i},'energywindow') || strcmpi(varargin{i},'energywindows')
        if nargin>i && isnumeric(varargin{i+1})
            energywindow = varargin{i+1} * 1e3;
            i=i+2;
        else
            error('problem with energywindow input argument')
        end
    elseif strcmpi(varargin{i},'depthdim') || strcmpi(varargin{i},'depth')
        if nargin>i && isnumeric(varargin{i+1})
            depthdim = varargin{i+1};
            i=i+2;
        elseif nargin>i && ischar(varargin{i+1}) && length(varargin{i+1})==1
            if strcmpi(varargin{i+1},'x') || strcmpi(varargin{i+1},'1')
                depthdim = 1;
                i=i+2;
            elseif strcmpi(varargin{i+1},'y') || strcmpi(varargin{i+1},'2')
                depthdim = 2;
                i=i+2;
            elseif strcmpi(varargin{i+1},'z') || strcmpi(varargin{i+1},'3')
                depthdim = 3;
                i=i+2;
            else
                error('problem with depthdim input argument')
            end
        else
            error('problem with depthdim input argument')
        end
    elseif strcmpi(varargin{i},'negdepth') || strcmpi(varargin{i},'negdepthflag')
        if nargin>i && (islogical(varargin{i+1}) || (isnumeric(varargin{i+1}) && (varargin{i+1}==1 || varargin{i+1}==0)))
            negdepthflag = logical(varargin{i+1});
            i=i+2;
        else
            negdepthflag = true;
            i=i+1;
        end
                
    else
        error(['Argument ID "',varargin{i},'" not recognized'])
    end
end
%default psft
if isempty(psft)
    tmp = load('psft.mat');
    psft = tmp.psft;
end

%shortcuts to column numbers
    indTrackID = 1;
    indParentID = 2;
    indStepNum = 3;
    indCharge = 4;
    indInitPos = 5:7;
    indFinalPos = 8:10;
    indTrackLen = 11;
    indStepLen = 12;
    indFinalE = 13;
    indInitE = 14;
    inddE = 15;
    
if isempty(who('CCDsize'))
    CCDsize = [726,1454];   %pixels
end
if variablesizeflag
    CCDsize = [4000,4000];  %pixels
end
pixsize = 10.5;         %microns

%% Mode: coincidence experiment simulation, OR single CCD.
%coordinates: center of CCD backplane is at (0,0,0); -z is into the detector.
if strcmp(inputmode,'coinc') || strcmp(inputmode,'single')
    %difference: coinc defaults to CCD depth in -z direction
    %            single defaults to CCD depth in +x direction
    
    %first, check input arguments.
    
    %set indDepth. This is a number 1, 2 or 3 corresponding to x, y, or z.
    % indDepth is the dimension along the drift of the charge carriers, perpendicular to the pixel plane.
    if ~isempty(who('depthdim'))
        %specified in input arguments
        indDepth = depthdim;
    else
        if strcmp(inputmode,'coinc')
            %coinc mode defaults to z
            indDepth = 3;
        else
            %single mode defaults to x
            indDepth = 1;
        end
    end
    indLateral = mod([indDepth+1,indDepth+2],3);
    indLateral(indLateral==0) = 3;
    %indLateral is whichver indices are not indDepth.
    % if indDepth = 3, then indLateral = [1,2];
    % if indDepth = 1, then indLateral = [2,3].

    if isempty(who('negdepthflag'))
        %if it is *not* specified in input arguments
        if strcmp(inputmode,'coinc')
            %coinc mode defaults to -
            negdepthflag = true;
        else
            %single mode defaults to +
            negdepthflag = false;
        end
    end
    
    %be explicit with what the max and min values of depth are.
    if negdepthflag
        mindepth = -.650;
        maxdepth = 0;
    else
        mindepth = 0;
        maxdepth = 0.650;
    end
    
    %zdepthflag is redundant with indDepth. But we use zdepthflag because DiffuseTrack3 code expects it.
    %If we ever need to use indDepth = 2, then we need to recode DiffuseTrack3.
    zdepthflag = (indDepth==3); %if indDepth==3, then zdepthflag = true; otherwise zdepthflag = false.
    
    if zcenteredflag        %Amy says this will never happen again....
        maxdepth = 0.325;
        mindepth = -0.325;
        %also, will adjust to be -0.65 < z < 0, for the diffusion code...
    end
    
    %separate the source photons
    f1 = find(M(:,indTrackID)==1 & M(:,indParentID)==0 & M(:,indStepNum)==1 & M(:,indCharge)==0);
    f1(end+1) = size(M,1)+1;    %for ease of indexing
    
    %progressbar increments by source photon.
    if progressbarflag
        progressbar(0);
    end
    
    %assign to output
    if length(f1)>2
        %multi-photon matrix: output will be a cell array of structures.
        %pre-allocate the cell array.
        out = cell(1,length(f1)-1);
%     else
%         %single photon matrix: output will be a structure. no preallocation needed
    end
    
    %now loop through source photons.
    for i=1:length(f1)-1    %each source photon
        Mtemp = M(f1(i):f1(i+1)-1,:);   %segment of matrix corresponding to this source photon
        
        %check event energy
        if Mtemp(find(Mtemp(:,indCharge),1),indInitE) < energywindow(1) || ...
                Mtemp(find(Mtemp(:,indCharge),1),indInitE) > energywindow(2)
            %out of energy window
            continue
        end
        
        if oldcodeflag
            D = DiffuseTrack2b(Mtemp,psft,650,false);   %noiseflag in DiffuseTrack2b should be hardcoded to false
            cheat = D.cheat;
            if isempty(noiselist)
                [T,E] = CCDsegment3(D.img,1e-4);
                if length(f1)==2
                    %single photon matrix: output a structure.
                    out.T = T;
                    out.E = E;
                else
                    %multi-photon matrix: output a cell array of structures.
                    out{i}.T = T;
                    out{i}.E = E;
                end
            else
                %multiple noise levels to simulate
                %must do everything multiple times
                for j=1:length(noiselist)
                    %add noise
                    imgtmp = AddCCDImageNoise(D.img,'blsig',noiselist(j));
                    %segment
                    %use variable threshold: 0.55 and blsig, added in quadrature
                    th = sqrt(0.55^2 + noiselist(j).^2);
                    [T,E] = CCDsegment3(imgtmp,th);
                    if length(f1)==2
                        %single photon matrix: one layer of cells, representing noise values
                        out{j}.T = T;
                        out{j}.E = E;
                        out{j}.cheat = cheat;
                        out{j}.blsig = noiselist(j);
                        out{j}.th = th;
                    else
                        %multi-photon matrix: two layers of cells
                        out{i}{j}.T = T;
                        out{i}{j}.E = E;
                        out{i}{j}.cheat = cheat;
                        out{i}{j}.blsig = noiselist(j);
                        out{i}{j}.th = th;
                    end
                end
            end
            continue
        end
        
        %initialize cheat structure.
        cheat.Edep = 0;     %so we can add each piece of energy deposition on, as it comes.
        
        %define image bounds.
        if variablesizeflag
            %use extents of the geant4 track, with buffer around
            buff = 20*pixsize*1e-3;  %20 pixels on each side. in mm.
            
            %ignore photons
            lgelectrons = Mtemp(:,indCharge)==-1;
            %ignore escapes
            lgcontained = Mtemp(:,indInitPos(indDepth))>mindepth & ...
                Mtemp(:,indInitPos(indDepth))<maxdepth;
            
            lg = lgelectrons & lgcontained;
            if ~any(lg)
                %no valid energy deposition here
                out = [];
                continue
            end
            
            %in g4 units (mm).
            minx = min([Mtemp(lg,indInitPos(indLateral(1))); Mtemp(lg,indFinalPos(indLateral(1)))]);
            maxx = max([Mtemp(lg,indInitPos(indLateral(1))); Mtemp(lg,indFinalPos(indLateral(1)))]);
            miny = min([Mtemp(lg,indInitPos(indLateral(2))); Mtemp(lg,indFinalPos(indLateral(2)))]);
            maxy = max([Mtemp(lg,indInitPos(indLateral(2))); Mtemp(lg,indFinalPos(indLateral(2)))]);
            
            %imagebounds = [minx, maxx; miny, maxy]
            imagebounds = [minx - buff, maxx + buff; miny - buff, maxy + buff]; %in mm
            imageboundspix = round(imagebounds *1e3./ pixsize);         %in pixels
            
            %redefine CCD size.
            CCDsize = (imageboundspix(:,2)-imageboundspix(:,1)+[1;1])';
            %CCDsize = [x,y]
        else
            %use real extents of CCD.
            %assume origin is at center of CCD in lateral dimensions.
            imagebounds = 1e-3*[-CCDsize(1)/2*pixsize, CCDsize(1)/2*pixsize; -CCDsize(2)/2*pixsize, CCDsize(2)/2*pixsize]; %in mm
            imageboundspix = round(imagebounds *1e3./ pixsize);  %in pixels
        end
        
        %adjust z if zcenteredflag. so that diffusetrack finds it in the right place.
        if zcenteredflag
            Mtemp(:,indInitPos(indDepth)) = Mtemp(:,indInitPos(indDepth)) - 0.325;
            Mtemp(:,indFinalPos(indDepth)) = Mtemp(:,indFinalPos(indDepth)) - 0.325;
            
            mindepth = mindepth - 0.325;
            maxdepth = maxdepth - 0.325;
        end
        
        %define image.
        img = zeros(CCDsize);
        
        %photon energy deposition (binding energy of electron)
        particlelist = unique(Mtemp(Mtemp(:,indCharge)==0,indTrackID));   %unique photon particleID
        for j=1:length(particlelist)
            %look for energy deposition inside CCD volume
            %do photon interactions one at a time, not all together for one photon (multiple compton scatter will be spread apart)
            lgCCDvol = Mtemp(:,indTrackID)==particlelist(j) & ...       current particle
                Mtemp(:,indFinalPos(indDepth)) > mindepth & ...         zmin
                Mtemp(:,indFinalPos(indDepth)) < maxdepth & ...         zmax
                Mtemp(:,indFinalPos(indLateral(1))) > imagebounds(1,1) & ...    xmin
                Mtemp(:,indFinalPos(indLateral(1))) < imagebounds(1,2) & ...    xmax
                Mtemp(:,indFinalPos(indLateral(2))) > imagebounds(2,1) & ...    ymin
                Mtemp(:,indFinalPos(indLateral(2))) < imagebounds(2,2) & ...    ymax
                Mtemp(:,inddE)>0;                                   % deposited energy
            fCCDvol = find(lgCCDvol);
            
            for k=1:length(fCCDvol)
                %diffuse this photon's energy deposition
                D = DiffuseTrack5(psft,'geant4',Mtemp(fCCDvol(k),:),'zdepth',zdepthflag,'negdepth',negdepthflag);
                %////debug
%                 disp(['Photon ',num2str(photonlist(j)),': Edep = ',num2str(D.cheat.Edep)])
                
                img = AddToImage(img,imageboundspix,D);
                
                %include photon energy deposition in total track energy deposition
                cheat.Edep = cheat.Edep + D.cheat.Edep;
            end
        end
        
        %electron energy deposition (multiple scattering, etc.)
        particlelist = unique(Mtemp(Mtemp(:,indCharge)==-1,indTrackID));    %unique electron particleID
        for j=1:length(particlelist)
            %look for energy deposition inside CCD volume (other deposition shouldn't be tracked, anyway)
            %  equal sign keeps steps that end at boundary.
            lgCCDvol = Mtemp(:,indTrackID)==particlelist(j) & ...       current particle
                Mtemp(:,indFinalPos(indDepth)) >= mindepth & ...         zmin
                Mtemp(:,indFinalPos(indDepth)) <= maxdepth & ...         zmax
                Mtemp(:,indFinalPos(indLateral(1))) >= imagebounds(1,1) & ...    xmin
                Mtemp(:,indFinalPos(indLateral(1))) <= imagebounds(1,2) & ...    xmax
                Mtemp(:,indFinalPos(indLateral(2))) >= imagebounds(2,1) & ...    ymin
                Mtemp(:,indFinalPos(indLateral(2))) <= imagebounds(2,2) & ...    ymax
                Mtemp(:,indInitPos(indDepth)) >= mindepth & ...         zmin
                Mtemp(:,indInitPos(indDepth)) <= maxdepth & ...         zmax
                Mtemp(:,indInitPos(indLateral(1))) >= imagebounds(1,1) & ...    xmin
                Mtemp(:,indInitPos(indLateral(1))) <= imagebounds(1,2) & ...    xmax
                Mtemp(:,indInitPos(indLateral(2))) >= imagebounds(2,1) & ...    ymin
                Mtemp(:,indInitPos(indLateral(2))) <= imagebounds(2,2) & ...    ymax
                Mtemp(:,inddE)>0;                                   % deposited energy
            
            if ~any(lgCCDvol)
                continue
            end
            %Mtemp2 is the portion of the matrix for this particleID and within the CCD active volume.
            Mtemp2 = Mtemp(lgCCDvol,:);
            %check for large gaps in between electron segments: if electron escapes CCD and then re-enters
            dx = Mtemp2(2:end,indInitPos)-Mtemp2(1:end-1,indInitPos);   %3-element vector [dx,dy,dz]
            dx2 = dx(:,1).^2 + dx(:,2).^2 + dx(:,3).^2;     %3D distance = sqrt(dx^2+dy^2+dz^2)
            dxthresh2 = (0.1).^2;      %mm^2                %hardcoded threshold here
            electronsplit = find(dx2>dxthresh2);        %find row indices where dx2 > dthresh2
            %make separate matrices for each segment
            %Mtemp3 is the matrix for one segment of an electron
            Mtemp3 = cell(1,1+length(electronsplit));
            if ~isempty(electronsplit)
                Mtemp3{1} = Mtemp2(1:electronsplit(1),:);
                for k=2:length(electronsplit)
                    Mtemp3{k} = Mtemp2(electronsplit(k-1)+1:electronsplit(k),:);
                end
                Mtemp3{end} = Mtemp2(electronsplit(end)+1:end,:);
            else
                %no splitting
                Mtemp3{1} = Mtemp2;
            end
            
            for k=1:length(Mtemp3)  %each segment of particle deposition
                %diffuse this electron's energy deposition
                D = DiffuseTrack5(psft,'geant4',Mtemp3{k},'zdepth',zdepthflag,'negdepth',negdepthflag);
                %////debug
%                     disp(['Electron ',num2str(particlelist(j)),', x=',num2str(Mtemp3{k}(1,5)),', segment #',num2str(k),': Edep = ',num2str(D.cheat.Edep)])
%                     disp(['Electron ',num2str(particlelist(j)),...
%                         ': Edep_mat = ',num2str(1e-3*sum(Mtemp(Mtemp(:,1)==particlelist(j),15))),...
%                         ' Edep_D = ',num2str(D.cheat.Edep)])
                
                img = AddToImage(img,imageboundspix,D);
                
                %if primary electron:
                if j==1
                    %add existing Edep to the cheat structure
                    Edep = cheat.Edep;
                    cheat = D.cheat;
                    cheat.Edep = cheat.Edep + Edep;
                else
                    %add Edep
                    cheat.Edep = cheat.Edep + D.cheat.Edep;
                end
            end
        end
        
        if isempty(noiselist)
            if ~nonoiseflag
                %add statistical and dark current noise
                img = AddCCDImageNoise(img);
                %segment
                [T,E] = CCDsegment3(img,0.55);
            else
                %no noise
                %save subset of image, using CCDsegment3
                [T,E] = CCDsegment3(img,1e-4);
            end
            
            
            %assign to output
            if length(f1)==2
                %single photon matrix: output a structure.
                out.T = T;
                out.E = E;
                if fullimageflag    %default false; can be set by input argument
                    out.img = img;
                end
                out.cheat = cheat;
            else
                %multi-photon matrix: output a cell array of structures.
                out{i}.T = T;
                out{i}.E = E;
                if fullimageflag
                    out{i}.img = img;
                end
                out{i}.cheat = cheat;
            end
        else
            %must do everything multiple times
            for j=1:length(noiselist)
                %add noise
                imgtmp = AddCCDImageNoise(img,'blsig',noiselist(j));
                %segment
                [T,E] = CCDsegment3(imgtmp,0.55);
                if length(f1)==2
                    %single photon matrix: one layer of cells, representing noise values
                    out{j}.T = T;
                    out{j}.E = E;
                    out{j}.cheat = cheat;
                    out{j}.blsig = noiselist(j);
                else
                    %multi-photon matrix: two layers of cells
                    out{i}{j}.T = T;
                    out{i}{j}.E = E;
                    out{i}{j}.cheat = cheat;
                    out{i}{j}.blsig = noiselist(j);
                end
            end
        end
        
        if progressbarflag
            %update progressbar
            progressbar(i/(length(f1)-1));
        end
    end
    
    if progressbarflag
        %make sure progressbar is closed
        progressbar(1);
    end
    
end


function img = AddToImage(img,imageboundspix,D)
    %truncate at edge of CCD active area
    trunc(1) = max(0, imageboundspix(1,1) - D.offsets(1));  %left
    trunc(2) = max(0,-imageboundspix(1,2) + D.offsets(1) + size(D.img,1)); %right
    trunc(3) = max(0, imageboundspix(2,1) - D.offsets(2));  %bottom
    trunc(4) = max(0,-imageboundspix(2,2) + D.offsets(2) + size(D.img,2)); %top
    
    x1 = -imageboundspix(1,1) + D.offsets(1) + (trunc(1)+1:(size(D.img,1)-trunc(2)));
    y1 = -imageboundspix(2,1) + D.offsets(2) + (trunc(3)+1:(size(D.img,2)-trunc(4)));
    %add to main image
    img(x1, y1) = img(x1,y1) + D.img(trunc(1)+1:end-trunc(2), trunc(3)+1:end-trunc(4));