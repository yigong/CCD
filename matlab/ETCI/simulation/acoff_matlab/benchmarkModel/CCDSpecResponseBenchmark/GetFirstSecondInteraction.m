function Rp = GetFirstSecondInteraction(R)

Rp = [];

F = find(R(:,2) == 0 & R(:,4)==0);

if F(1) == 1 & length(F) > 2
    Rp = R(F(2:3),:);
    
end

 