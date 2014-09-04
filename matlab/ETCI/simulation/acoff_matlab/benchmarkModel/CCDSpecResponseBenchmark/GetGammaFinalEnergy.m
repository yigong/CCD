function [Egam] = GetGammaFinalEnergy(R) 

% Processes the Eneryg Spectrum of the Ge detector

Re = [];
Egam = 0;

F = find(R(:,2) == 0 & R(:,4)==0);

if F(1) == 1 & length(F) > 2
    Re = R(F(2:length(F)),:);
end

E = 0;

if isempty(Re)
    return;
elseif size(Re,1) == 1
    E = Re(:,14)-Re(:,13);
    return;
end
%Calculate Track Energy in Ge
E1 = Re(:,14); E2 = Re(:,13);
dE = E2(1:end-1)-E2(2:end);
dE = [E1(1)-E2(1);dE];
E = sum(dE);

% Temp Output
Egam = E;
