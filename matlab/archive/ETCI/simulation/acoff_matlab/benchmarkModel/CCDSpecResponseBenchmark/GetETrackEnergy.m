function [E,dE] = GetETrackEnergy(Re)

E = 0;

if isempty(Re)
    return;
elseif size(Re,1) == 1
    E = Re(:,14)-Re(:,13);
    return;
end

%Calculate Track Energy
E1 = Re(:,14); E2 = Re(:,13);
dE = E2(1:end-1)-E2(2:end);
dE = [E1(1)-E2(1);dE];
E = sum(dE);