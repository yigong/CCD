function bool = QualifyGammaTrack(R)

F = find(R(:,2) == 0 & R(:,4) == 0);

if sum(R(F(1:end-1),13)-R(F(2:end),14)) == 0
    bool = true;
    return;
else
    bool = false;
    return;
end
    