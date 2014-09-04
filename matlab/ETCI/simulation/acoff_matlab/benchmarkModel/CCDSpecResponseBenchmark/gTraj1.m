function out = gTraj1(S)

%%%%%%%%%%%%NEW New Matrix%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Column 1: Track ID
% Column 2: Step Number
% Column 3: Parent ID
% Column 4: Charge
% Column 5,6,7: Initial position (x,y,z)
% Column 8,9,10: Final position (x,y,z)
% Column 11: Track Length
% Column 12: Step Length
% Column 13: Final Energy (E)
% Column 14: dE
stepNum = S(:,1);
ld=length(stepNum);
f = [];
% Finds if there are any step ids other than 2 tracks start:
% Returns index of each starting electron track.
for i=1:ld
    %if abs(dx(i))>.01 && abs(dy(i))>.01 && abs(dz(i))>.01
    if stepNum(i)~=2
        f=[f; i+1];
    end
end

% Print out what the trackID is instead of 2
trackID = [];
for i=1:length(f)
   trackID = [trackID; S(f(i)-1)];
end
out = trackID;
