% Pt sources and their transforms: 
% From Experiment Reference Frame to Simulation reference frame (CCD depth
% in +x-axis

%
% pBeam
%

%
% Cs-137 Pt Source
%

%
% Ba-133 Pt Sources
%

% Angle-1 
% Start (-187.07 29.02 140.85) [mm]
% A
x1 = -187.07;
y1 = 29.02;
z1 = 140.85;

% CCD Front Face (0 0 0) [mm]
%B
x0 = 0; y0 = 0; z0 = 0.;
%C
xC =0; yC =1.5; zC=1.5; % [mm]

A = [-187.07 29.02 140.85];
B = [0 0 0];
C = [0 1.5 1.5];

aP = [B' C']%B-C
aV = dist(aP);
a = max(max(aV)) %[mm]  Distance between tri's A and C

bP = [A' C'];
b = max(max(dist(bP))) % [mm] Distance

cP = [A' B'];
c = max(max(dist(cP)))

ang_C = acosd((a^2+b^2-c^2)/(2*a*b));
ang_A = acosd((c^2+b^2-a^2)/(2*c*b))

x2 = (x0-x1);
y2 = (y0-y1);
z2 = (z0-z1);

aMag = sqrt(x2*x2+y2*y2+z2*z2);
ux = (x2/aMag)
uy = y2/aMag
uz = z2/aMag



% Openening Angle:


%
% Co-60 Pt Sources
%


disp(['Beam Source Start Location (' num2str(x1) ' ' num2str(y1) ' ' num2str(z1) ')'])
disp(['In the Unit Direction Vector of (' num2str(ux) ' ' num2str(uy) ' ' num2str(uz) ')'])
