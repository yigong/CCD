function bool = QualifyETrack(Re,minEnergy)

bool = false;

if isempty(Re)
    return;
end

if GetETrackEnergy(Re) < minEnergy
    return;
end

%Find escape tracks
x = Re(:,8); y = Re(:,9); z = Re(:,10);
%dx = abs(x(2:end)-x(1:end-1));
%if sum(dx>0.5) % Greater than 0.5mm
%    return;
%end

%Find Tracks not spatially local
dX = abs(max(x)-min(x));
dY = abs(max(y)-min(y));
dZ = abs(max(z)-min(z));

CT = 2; %2mm
if dX>CT || dY>CT || dZ>CT
    return;
end



bool = true;
    