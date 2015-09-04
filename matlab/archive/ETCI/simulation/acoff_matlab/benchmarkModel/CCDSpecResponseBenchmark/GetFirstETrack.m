function RTrk = GetFirstETrack(R)

RTrk = [];

F = find(R(:,2) == 0 & R(:,4)==0);

if F(1) == 1 & length(F) > 1
    F = F(2);
end


%Fe = find(R(:,2)==1 & R(:,4)==-1 & R(:,5)==R(F,8) & R(:,6)==R(F,9) & R(:,7)==R(F,10))
Fe = find(R(:,1)==2 & R(:,2)==1 & R(:,3)==1 & R(:,4)==-1);

if isempty(Fe)
    disp('Warning: First Electron Track Not Found');
    return;
end
Fe = Fe(1);

IDe = R(Fe,1);

RTrk = R(find(R(:,1) == R(Fe,1)),:);

