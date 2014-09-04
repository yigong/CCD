function [Event] = RunTracks7(Res,minEThresh, dEdx_ref, psft, depthCordinates)
%[Event] = RunTracks3(Res,minEThresh, dEdx_ref, psft, runNPS, depthCordinates)
% UPDATED: 06/04/12 - Coffer - Adding ARM and Gamma-ray tracking to Outputs
% load(filename);
% 
% OUTPUTS:
% Event{i} = struct('trackM',R2,'out',out, 'Etot', Etot,
% 'Edep',Edep,'Eesc', Eesc, 'multiplicity', multiplicity, 'geEngDep',{En},'DiffErr', D,'boolenQuant',[QET,QGT]);
% out:
%       out.
%       out.
% ADD:


if isempty(Res)
    return;
end

%load('../psft_650um.mat') % psft_650um.mat')         % 650um psft Table
%load('../dEdx_ref.mat')     % Silicon dEdx ref table
%tunc = load('TrackUnc.mat');

% Detector Construction Variables
th = 650;
lengthRes = length(Res);

% Ge Detecor Configuration
% Ge Depth Location (z < -40 )
initialDepthPos = 7; % z
minDetectorDepthPosition = -40;  % [mm]
maxDetectorDepthPosition = -120; % [mm]


% Declare Variables
 
En = [];
RDF = [];
RHT = [];
%QET = [];
%QGT = [];
Egam = [];
Etot = [];
Edep = [];
Eesc = [];
multiplicity = [];
tic

for i = 1:length(Res)
   
    if (mod(i,10) == 0)
        toc
        disp('*********************');
        disp(i)
        tic
    end
    % Track Matrix
    R2 = Res{i};
    %Re = GetFirstETrack(R2);  
    %Rp2 = GetFirstSecondInteraction(R2);
    %En = GetETrackEnergy(Re);
    
    %QET = QualifyETrack(Re,20000);
    %QGT = QualifyGammaTrack(R2);
    
    clear D;
    D.err = 'No Attempt';
  %  clear A;
  %  A.err = 'No Attempt';
   
        % Run Diffuse Track and Find and Save Electron Track (CCD-Detector)
     %   try
            %D = DiffuseTrack4(R2,psft,th,0);
            % Need this Etot and Edep for determining 'Escapes' later.
            % Records CCD Info
            [Etot,Edep,Eesc,multiplicity] = GetEnergiesFromGeant(R2);

            % For NEW Geant4TrackHandling: [2000 2000] number of CCD pixels
            out = Geant4TrackHandling7(R2,'coinc','psft',psft,'CCDsize',[2000,2000],'bypassPhotonErrorFlag', 1, 'driftdim', 'z', 'depthcoordinates', depthCordinates);
            % Successfuly Diffused
            out.err = [];
            D.err = [];
%         catch e
%             D.err = 'Error Diffusing Track';
%             D.errmsg = e.message;
%             D.errStack = e.stack;
%             D.errIdentifier = e.identifier;
%             D.errCause = e.cause;
%             out.err = e;
%         end
        
%    % Only Useful for Coinc Data Sets  
        % Find and Save Gamma-Ray Tracking Info (Ge-Detector)
%         try 
%             GeDetector = GetGeGammaEnergies(R2,initialDepthPos,minDetectorDepthPosition,maxDetectorDepthPosition);
%             GeDetector.err = [];
%         catch e
%             GeDetector.err = e;
%         end
%         
%         % Run and Save ARM Calculations
%         try
%             % ARM Calculation Using CCD and Ge Detector Info
%             armInfo = calculateARM6(out,GeDetector);
%             armInfo.err = [];
%         catch e
%             armInfo.err = e;
%         end
        
        
% Run Brian's Electron Trajectory Algorithm
%    % Might be useful to add an Electron 'ARM' based just on The Tracking
%    % Algorithm at some point... 
%    %OLD %if QualifyETrack(Re,minEThresh)  & ~isempty(Rp2) %& QualifyGammaTrack(R2)
%           
%    % OK --     
%         clear A;
%         A.err = 'Post-Diffuse No Attempt';
%         
%         if isempty(D.err)
%             try
%                 A = HybridTrack_1(D.img,0.08,dEdx_ref,[],0);
%                 A.err = [];
%             catch e
%                 A.err = 'Error in Hybrid Track';
%                 A.errmsg = e.message;
%             end
% 
%         end
%     %end
    
    % Save Variables
    %Event(i) = struct('trackM',R2,'diffIMG',{D},'hybridInfo',{A},...
    %    'geEngDep',{En},'boolenQuant',[QET,QGT,Egam]);
    Event{i} = struct('trackM',R2,'out',out, 'Etot', Etot, 'Edep',Edep, ...
        'Eesc', Eesc, 'multiplicity', multiplicity,  ...
        'DiffIMG', D, 'depthCordinates', depthCordinates);  %'geDetector',GeDetector,'armInfo', armInfo, ...'boolenQuant',[QET,QGT]);
    
    % Reset Variables
    En = []; RDF = []; RHT = [];  Egam = [];
    Etot = []; Edep = []; Eesc = []; multiplicity = []; out = [];
    %GeDetector = [];    % QET = []; QGT = []; % Old Quantities %
    clear R2 % Re Rp2 
end

toc

% F = strfind(filename,'.mat');
% sname = [filename(1:(F-1)),'_TRK.mat'];
% save(sname,'R','RDF','RHT','En','QET','QGT');
% RTrk=0;
        
