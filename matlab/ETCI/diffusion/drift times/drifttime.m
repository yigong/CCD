function drifttable = drifttime(step,Vdep)
%function drifttable = drifttime(step,Vdep)
%
% Edep is depletion voltage in V
% step is step size in um
%
% Returns table:
%
%Col 1: depth(z)    (0 = backside, 650 = pixels)    um
%Col 2: E(z)        E(z) = ezcoef * z               V/cm
%Col 3: Vd          Vd(E) = interp1(Vd140)          cm/s
%Col 4: Vd                                          um/s
%Col 5: Td          Td = 1/Vd * dz                  s
%
% Col 5 is differential drift time over dz
% 
% Then use totaldrifttime.m


T = 140;
th = 650;
Vapp = 200;





%Vd =  a*x/((1+(x/b)^c)^(1/c))
% with x representing electric field (V/cm)

% %from fit codes (fit110Kcustom.m):
% T110a = 4449.1120491440479;
% T110b = 2691.1604838621606;
% T110c = 0.76934879612066231;
% 
% T160a = 2274.1443149396855;
% T160b = 4579.4489122565683;
% T160c = 0.84964350356874552;

%from Roland's code (better fit):
T110a = 3850.3818;
T110b = 2776.3685466505467;
T110c = 0.90727407870716548;

T160a = 2056;
T160b = 4505.5247670943809;
T160c = 0.98909318158869919;


E = 0:0.1:7000;
Vd110 = T110a*E./((1+(E/T110b).^T110c).^(1/T110c));
Vd160 = T160a*E./((1+(E/T160b).^T160c).^(1/T160c));
% u110 = Vd110./E;
% u160 = Vd160./E;

%logarithmic interpolation
log110 = log(Vd110);
log160 = log(Vd160);
log140 = log110 + (T-110)./(160-110).*(log160 - log110);
Vd140 = exp(log140);
Vd140(1) = 0;

drifttable(:,1) = 0:step:th;

%Col 1: depth(z)    (0 = backside, 650 = pixels)    um
%Col 2: E(z)        E(z) = [see below]               V/cm
%Col 3: Vd          Vd(E) = interp1(Vd140)          cm/s
%Col 4: Vd                                          um/s
%Col 5: Td          Td = 1/Vd * dz                  s

%E(z) = -2/650^2 *Vdep*z - 1/650*(Vapp - Vdep) = V/um
drifttable(:,2) = (2/th^2*drifttable(:,1)*Vdep + (Vapp - Vdep)/th) * 10000;    %in V/cm now

drifttable(:,3) = interp1(E,Vd140,drifttable(:,2));
drifttable(:,4) = drifttable(:,3).*10000;   %cm to um
drifttable(:,5) = step ./ drifttable(:,4);
if drifttable(1,4)==0
    drifttable(1,5) = drifttable(2,5);  %linear extrapolation...
end


% 
% plot(E,u110)
% hold on
% plot(E,u160)
