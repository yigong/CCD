function errTest = calculateARM(electronTRKArray)
% errTest = calculateARM(electronTRKArray)

load(electronTRKArray);
n = length(QET);
errTest = {};

% Angle between the known direction of the initial and measure direction of
% the scattered gamma ray (by scatter geometry)
%Pencil Beam info: a-->b, ei - Direction of the inital gamma-ray
a = [-37.5168, 92.8575, 0];
b = [-28.4701, 70.46597, 0];
ba = (b-a);
magBA = sqrt(dot(ba,ba));
ei = ba./magBA;
eiMag = sqrt(dot(ei,ei));   % Test to make sure its Normalized

for i = 1:n
    
    if ~isempty(RDF{1,i}.err)   % Make sure the track diffued correctly
        errTest.err(i) = 0;       % Zero Array if failed
        return
    else
        errTest.err(i) = 1;
    end

    % True Theta
    thetaTrue = 90;    % [deg]
    
    % ee - Direction of the scattered electron 
    eTrack1 = GetFirstETrack(R{1,i});  
    ee1 = [eTrack1(1,5), eTrack1(1,6), eTrack1(1,7)];
    ee2 = [eTrack1(1,8), eTrack1(1,9), eTrack1(1,10)];
    ee21 = (ee2-ee1);
    magEE21 = sqrt(dot(ee21,ee21));
    ee = ee21./magEE21;
    
    % eg - Direction of the scattered gamma-ray
    gTrack1 = GetFirstSecondInteraction(R{1,i});  
    eg1 = [gTrack1(1,8), gTrack1(1,9), gTrack1(1,10)];
    eg2 = [gTrack1(2,5), gTrack1(2,6), gTrack1(2,7)];
    eg21 = (eg2-eg1);
    magEG21 = sqrt(dot(eg21,eg21));
    eg = eg21./magEG21;
    
    E1 = RDF{1,i}.cheat.Etot;  % Need to Fix this value to correct Etrack value in new Diffusion code
    Ei = 662.0;           % [keV]
    E2 = (Ei-E1);         % [keV]
    eMass = 510.998910;   % [keV]
    
    % Calc Cos(theta)
    cosThetaMeasured = (1- eMass./E2+eMass./(E2+E1));
    thetaMeasured = acos(cosThetaMeasured);
    thetaMeasuredDeg = (180/pi).*thetaMeasured;
    errTest.theta(i) = thetaMeasuredDeg;
    
    % Photon-ARM
    errTest.gARM(i) = acos(dot(ei,eg))-thetaMeasured;
    errTest.gARMdeg(i) = (180/pi).*(errTest.gARM(i));
    % Electron-ARM
    errTest.eARM(i) = acos(dot(ee,eg))-thetaMeasured;
    errTest.eARMdeg(i) = (180/pi).*(errTest.eARM(i));
    % Scatter Plane Deviation - SPD
    errTest.SPD(i) = acos(dot(cross(eg,ei),cross(eg,ee)));
    errTest.SPDdeg(i) = (180/pi).*errTest.SPD(i);
end

