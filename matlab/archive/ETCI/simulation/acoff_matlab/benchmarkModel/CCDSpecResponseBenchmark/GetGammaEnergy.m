function Egam = GetGammaEnergy(RTrk) 

% Processes the Eneryg Spectrum of the Ge detector

Rp = [];

F = find(R(:,2) == 0 & R(:,4)==0);

if F(1) == 1 & length(F) > 2
    Rp = R(F(2:length(F)),:);
end

