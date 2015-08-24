function ARMinfo = calculateARM5(Event)
% errTest = calculateARM(electronTRKArray)
% Event Data Structure: Event{1,i} from convertBinary5.m
% OUTPUTS:
% %CCD
%   ARMinfo.E_CCDimg
% %Germanium Detector
%   ARMinfo.E_GeDet_dE
%   ARMinfo.EtotalGe 
%   ARMinfo.E_CCD_Plus_Ge
% %Angles
%   ARMinfo.theta
% %ARM
% %Photon-ARM
%   ARMinfo.gARM
%   ARMinfo.gARMdeg
% %Electron-ARM
%   ARMinfo.eARM 
%   ARMinfo.eARMdeg 
% %Scatter Plane Deviation - SPD
%   ARMinfo.SPD
%   ARMinfo.SPDdeg

%load(electronTRKArray);
%n = length(Event);
ARMinfo = {};
ARMinfo.err = {};
n =1; 

% Initalize with Empty Arrays:
ARMinfo.errBool = [];
ARMinfo.E_CCDimg = [];
ARMinfo.E_CCD_Plus_Ge = [];
ARMinfo.E_GeDet_dE = [];
ARMinfo.EtotalGeTemp = [];

% Purly ARM Output Set
ARMinfo.theta = [];
ARMinfo.gARM = [];
ARMinfo.gARMdeg = [];
ARMinfo.EtotalGe = [];
ARMinfo.E_ARM_CCD = [];
%ARMinfo.EtotalGe = 0;
%ARMinfo.E_ARM_CCD = 0;
ARMinfo.eARM = [];
ARMinfo.eARMdeg = [];
ARMinfo.SPD = [];
ARMinfo.SPDdeg = [];



% Angle between the known direction of the initial and measure direction of
% the scattered gamma ray (by scatter geometry)
%Pencil Beam info: a-->b, ei - Direction of the inital gamma-ray
% October 2011 Simulations:
% x1 = -925.2782*mm,     //-1190.6*mm,   //8.094*mm,
% y1 = -8.3933*mm,       //417.6*mm,     //4.266*mm,
% z1 = 374.5824*mm; 
% x0 = 1.906*mm,
% y0 = -4.221*mm,
% z0 = 0.;

% Beam Vector Calculation
% a = [-37.5168, 92.8575, 0];
a = [ -925.2782, -8.3933, 374.5824];
% b = [-28.4701, 70.46597, 0];
b = [ 1.906, -4.221, 0];

ba = (b-a);
magBA = sqrt(dot(ba,ba));
ei = ba./magBA;
%eiMag = sqrt(dot(ei,ei));   % Test to make sure its Normalized

    try 
        % Find Scattered Electron Info
        % ee - Direction of the scattered electron 
        % eTrack1 = GetFirstETrack(R{1,i});  
        %eTrack1 = GetFirstETrack(Event{1,i}.trackM); 
        %ee1 = [eTrack1(1,5), eTrack1(1,6), eTrack1(1,7)];
        %ee2 = [eTrack1(1,8), eTrack1(1,9), eTrack1(1,10)];
        ee1 = [Event.out.cheat.x0(1,:)];
        ee2 = [Event.out.cheat.x(2,:)];
        ee21 = (ee2-ee1);
        magEE21 = sqrt(dot(ee21,ee21));
        ee = ee21./magEE21;
        % Find Scattered Gamma-Ray Info
        % eg - Direction of the scattered gamma-ray
        %gTrack1 = GetFirstSecondInteraction(Event{1,i}.trackM);  
        %eg1 = [gTrack1(1,8), gTrack1(1,9), gTrack1(1,10)];
        %eg2 = [gTrack1(2,5), gTrack1(2,6), gTrack1(2,7)];
        % Ge-Detector Position and Energies:
        geDetector = GetGeGammaEnergies(Event.trackM);       
        
        %Calculate Gamma-Ray Scattering Parameters
        % Gamma-ray Direction Vector starts with Electron Track Start Position
        eg1 = [Event.out.cheat.x0]; % 
        eg2 = geDetector.gamGePosition(1,:);
        eg21 = (eg2-eg1);
        magEG21 = sqrt(dot(eg21,eg21));
        eg = eg21./magEG21;
        

        
        %ARMinfo.E_CCD_Plus_Ge = ARMinfo.E_CCDimg + geDetector.EtotalGe; % [keV]
        ARMinfo.E_CCDimg = Event.out.E;
        ARMinfo.E_GeDet_dE = geDetector.dE;
        ARMinfo.EtotalGeTemp = geDetector.EtotalGe;
    
    % Run ARM Calculations
        % Calc Cos(theta)
        %cosThetaMeasured = (1- eMass./E2+eMass./(E2+E1));
        % Only Calculate ARM info if All Data is Avaiable:
        if (geDetector.err==1) 
            disp('Bad Ge-Det Info')
            return
        elseif isempty(ARMinfo.E_CCDimg) 
            disp('Bad CCD-Det Info')
            return
            
        elseif isempty(geDetector.EtotalGe)
            disp('Bad Ge-Det 2 Info')
            return
        else
            ARMinfo.EtotalGe = geDetector.EtotalGe;
            ARMinfo.E_ARM_CCD = (Event.out.E);
            % Energies
            E1 =(Event.out.E);          % [keV]
            E2 = ARMinfo.EtotalGeTemp;  % [keV]
            eMass = 510.998910;         % [keV]
            
            ARMinfo.E_CCD_Plus_Ge = ARMinfo.EtotalGe + ARMinfo.E_ARM_CCD; % [keV]
            cosThetaMeasured = 1 - (eMass*E1)./(E2.*(E1+E2));
            thetaMeasured = acos(cosThetaMeasured);
            thetaMeasuredDeg = (180/pi).*thetaMeasured;
            ARMinfo.theta = thetaMeasuredDeg;

            % Photon-ARM
            %mu(i,1) = dot(beamvec,dpos) ./ (norm(beamvec) * norm(dpos));
            ARMinfo.positionTheta = acos(dot(ei,eg)); 
            ARMinfo.gARM = ARMinfo.positionTheta-thetaMeasured;
            ARMinfo.gARMdeg = (180/pi).*(ARMinfo.gARM);
            
            % Electron-ARM
            ARMinfo.eARM = acos(dot(ee,eg))-thetaMeasured;
            ARMinfo.eARMdeg = (180/pi).*(ARMinfo.eARM);
            
            % Scatter Plane Deviation - SPD
            ARMinfo.SPD = acos(dot(cross(eg,ei),cross(eg,ee)));
            ARMinfo.SPDdeg = (180/pi).*ARMinfo.SPD;
 
            ARMinfo.errBool = [ARMinfo.errBool 0];
        end
        
    catch e
        ARMinfo.err{n} = e;
        n = n+1;
        ARMinfo.errBool = [ARMinfo.errBool 1];
    end
        
