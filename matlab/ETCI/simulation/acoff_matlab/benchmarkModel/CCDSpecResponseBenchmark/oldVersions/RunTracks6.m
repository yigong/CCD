function [Event] = RunTracks6(Res,minEThresh, dEdx_ref, psft)
%[Event] = RunTracks3(Res,minEThresh, dEdx_ref, psft)
%load(filename);

if isempty(Res)
    return;
end

%load('../psft_650um.mat') % psft_650um.mat')         % 650um psft Table
%load('../dEdx_ref.mat')     % Silicon dEdx ref table
%tunc = load('TrackUnc.mat');

th = 650;
lengthRes = length(Res)


% Declare Variables
 
En = [];
RDF = [];
RHT = [];
En = [];
QET = [];
QGT = [];
Egam = [];
Etot = [];
Edep = [];
Eesc = [];
multiplicity = [];

for i = 1:length(Res)
   
    disp('*********************');
    disp(i)
    R2 = Res{i};
    Re = GetFirstETrack(R2);  
    Rp2 = GetFirstSecondInteraction(R2);
    En = GetETrackEnergy(Re);
    QET = QualifyETrack(Re,20000);
    QGT = QualifyGammaTrack(R2);
    
    % Determine Gamma Ray energy deposited in Ge
    %   if QGT(i)

%Egam(i) = GetGammaFinalEnergy(R2);
  %  else
        %Egam(i) = 0;
  %  end
        
    clear D;
    D.err = 'No Attempt';
    clear A;
    A.err = 'No Attempt';
    
    %if QualifyETrack(Re,minEThresh)  & ~isempty(Rp2) %& QualifyGammaTrack(R2)
   
        %clear D
        
        try
            %D = DiffuseTrack4(R2,psft,th,0);
            [Etot,Edep,Eesc,multiplicity] = GetEnergiesFromGeant(R2);
            %out{i} = Geant4TrackHandling(R2,'coinc','psft',psft);
            %For Coincident experiment geometry set-up
            % OLD Geant4TrackHandling:
            %out = Geant4TrackHandling(R2,'coinc','psft',psft);
            %For Brian's single CCD e_Tracks geant4 
            %out = Geant4TrackHandling(R2,'single','psft',psft);
            % For NEW Geant4TrackHandling: [2000 2000] number of CCD pixels
            out = Geant4TrackHandling5(R2,'coinc','negdepth',true,'CCDsize',[2000,2000],'psft',psft);
            
            if ~isempty(out) && ~isempty(out.E); 
                dE(i) = out.E(1) - Event{i}.Edep(1); 
            else
                dE(i) = nan;
            end
            %D = out.
            
            D.err = [];
        catch e
            D.err = 'Error Diffusing Track';
            D.errmsg = e.message;
            D.errStack = e.stack;
            D.errIdentifier = e.identifier;
            D.errCause = e.cause;
        end
        
       
%     if QualifyETrack(Re,minEThresh)  & ~isempty(Rp2) %& QualifyGammaTrack(R2)
%           
%         
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
%     end
    
    % Save Variables
    %Event(i) = struct('trackM',R2,'diffIMG',{D},'hybridInfo',{A},...
    %    'geEngDep',{En},'boolenQuant',[QET,QGT,Egam]);
    Event{i} = struct('trackM',R2,'out',out, 'Etot', Etot, 'Edep',Edep,'Eesc', Eesc, 'multiplicity', multiplicity, 'geEngDep',{En},'boolenQuant',[QET,QGT]);
    
    % Reset Variables
    En = []; RDF = []; RHT = []; En = [];QET = []; QGT = []; Egam = [];
    Etot = []; Edep = []; Eesc = []; multiplicity = []; out = [];
    clear R2 Re Rp2
end
% F = strfind(filename,'.mat');
% sname = [filename(1:(F-1)),'_TRK.mat'];
% save(sname,'R','RDF','RHT','En','QET','QGT');
% RTrk=0;
        
