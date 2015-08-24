% Pencil Beam Parameters and Calculation for Coincident Experiment 2010
%

x1 = -925.2782;
y1 = -8.3933;
z1 = 374.5824;

x0 = 1.906;
y0 = -4.221;
z0 = 0.;

x2 = (x0-x1);
y2 = (y0-y1);
z2 = (z0-z1);

aMag = sqrt(x2*x2+y2*y2+z2*z2);
ux = (x2/aMag);
uy = y2/aMag;
uz = z2/aMag;

disp(['Beam Source Start Location (' num2str(x1) ' ' num2str(y1) ' ' num2str(z1) ')'])
disp(['In the Unit Direction Vector of (' num2str(ux) ' ' num2str(uy) ' ' num2str(uz) ')'])
