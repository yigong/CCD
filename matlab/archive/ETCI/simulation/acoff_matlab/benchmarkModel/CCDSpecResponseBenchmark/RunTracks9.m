function [Event] = RunTracks9(Res,minEThresh, dEdx_ref, psft, depthCordinates,radDecayFlag)
%[Event] = RunTracks3(Res,minEThresh, dEdx_ref, psft, runNPS, depthCordinates)
% UPDATE: 06/04/12 - Coffer - Adding ARM and Gamma-ray tracking to Outputs
% UPDATE TWO: 11/27/2013 - Coffer - Modifying for the New Tracking Matrix
% UPDATE 3: 2/7/2014 - Adding Error catching aroung G4TrackingHandling. (out.err catches error info)
% load(filename);
% 
% OUTPUTS:
% Event{i} = struct('trackM',R2,'out',out, 'Etot', Etot,
% 'Edep',Edep,'Eesc', Eesc, 'multiplicity', multiplicity, 'geEngDep',{En},'DiffErr', D,'boolenQuant',[QET,QGT]);
% out: (still in flux)
%       out.
%       out.
% 

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
%initialDepthPos = 7; % z
initialDepthPos = 5;  % x
minDetectorDepthPosition = 0;  % [mm]
maxDetectorDepthPosition = 0.650; % [mm]


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
    
    %clear out;
    clear err1
    err1 = 'No Attempt';
  %  clear A;
  %  A.err = 'No Attempt';
   
        % Run Diffuse Track and Find and Save Electron Track (CCD-Detector)
        try
            %D = DiffuseTrack4(R2,psft,th,0);
            % Need this Etot and Edep for determining 'Escapes' later.
            % Records CCD Info
            % Obsolete: 
            %[Etot,Edep,Eesc,multiplicity] = GetEnergiesFromGeant8(R2);

            % For NEW Geant4TrackHandling: [2000 2000] number of CCD pixels
            %out = Geant4TrackHandling10(R2,'coinc','psft',psft,'CCDsize',[2000,2000],'bypassPhotonErrorFlag', 1, 'driftdim', 'x', 'depthcoordinates', depthCordinates,'RadDecay', radDecayFlag);
            out = Geant4TrackHandling10(R2,'single','psft',psft,'variablesize',true,'bypassPhotonErrorFlag', true, 'driftdim', 'x', 'depthcoordinates', depthCordinates,'RadDecay', radDecayFlag);
            %variableSizeFlag
            % Successfuly Diffused
            err1 = [];
            %D.err = []; 
         catch e
%             D.err = 'Error Diffusing Track';
%             D.errmsg = e.message;
%             D.errStack = e.stack;
%             D.errIdentifier = e.identifier;
%             D.errCause = e.cause;
             err1 = 'Failed to Diffuse properly';   
             out.err = e;
         end
        
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
%     Event{i} = struct('trackM',R2,'out',out, 'Etot', Etot, 'Edep',Edep, ...
%         'Eesc', Eesc, 'multiplicity', multiplicity,  ...
%         'DiffIMG', D, 'depthCordinates', depthCordinates);  %'geDetector',GeDetector,'armInfo', armInfo, ...'boolenQuant',[QET,QGT]);
    %trackM{i} =  struct('trackM',R2);
    Event{i} = struct('out',out, 'err1',err1,'depthCordinates', depthCordinates);    
    % Reset Variables
    En = []; RDF = []; RHT = [];  Egam = [];
    Etot = []; Edep = []; Eesc = []; multiplicity = []; out = [];
    %GeDetector = [];    % QET = []; QGT = []; % Old Quantities %
    clear R2 out % Re Rp2 
end

% Display final time:
disp('*********************');
disp(i)
toc

% F = strfind(filename,'.mat');
% sname = [filename(1:(F-1)),'_TRK.mat'];
% save(sname,'R','RDF','RHT','En','QET','QGT');
% RTrk=0;
        
