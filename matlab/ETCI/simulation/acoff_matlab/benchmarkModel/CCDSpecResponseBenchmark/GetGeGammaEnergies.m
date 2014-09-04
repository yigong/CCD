function GeDetector = GetGeGammaEnergies(RTrk) 
% GeDetector = GetGeGammaEnergies(RTrk) 
% RTrk = Event{1,i}.trackM   : Event's Tracking Matrix 

% Record the Gamma Energies in the Ge Detector 
%
% Ge Detecor Configuration
% Ge Depth Location (z < -40 )
initialDepthPos = 7; 
minDetectorDepthPosition = -40;  % [mm]
maxDetectorDepthPosition = -120; % [mm]

GeDetector = {};

Re = [];
Egam = 0;
% Outputs
GeDetector.dE = [];             %   [keV]
GeDetector.EtotalGe = [];       %   [keV]
GeDetector.gamGeTrack = [];     %   [mm]
GeDetector.gamGePosition = [];  %   [mm]
EtotalGe = [];


% Saves the Gamma-Ray  Track from Event{}.trackM
F = find(RTrk(:,2) == 0 & RTrk(:,4)==0);
if F(1) == 1 & length(F) > 2
    Re = RTrk(F(1:length(F)),:);
end


if isempty(Re)
    return;
end

% Find the Gamma-Ray Tracking Info through the Ge-Detector
F2 = find(Re(:,2) == 0 & Re(:,4)==0 & ...
     (Re(:,initialDepthPos)<minDetectorDepthPosition) & ...
      (Re(:,initialDepthPos)>maxDetectorDepthPosition) );

Re2 = Re(F2(1:length(F2)),:);


% Find the Gamma-Ray Interaction Locations in Ge-Detector
F3 = find(abs(Re2(:,13)-Re2(:,14))>0);
    Re3 = Re2(F3(1:length(F3)),:);
 
    
    %Calculate Track Energy in Ge
    E1 = Re3(:,14); 
    E2 = Re3(:,13);
    Edep = Re3(:,15);
    gamPosition = [Re3(:,8) Re3(:,9) Re3(:,10)];
    dE = (E1(1:end)-E2(1:end)-Edep(1:end));
    EtotalGe = sum(dE);
    
if isempty(EtotalGe)
    GeDetector.err =1;
    return
else
    % Outputs
    %GeDetector.Edep = Edep;
    GeDetector.err =0;
    GeDetector.dE = dE./(1e3);              %   [keV]
    GeDetector.EtotalGe = EtotalGe./(1e3);   %   [keV]
    GeDetector.gamGeTrack = Re3;            %   [mm]
    GeDetector.gamGePosition = gamPosition; %   [mm]
end

    

%GeDetector.E_Ge1 = ;
%GeDetector.E_Getot = ;
