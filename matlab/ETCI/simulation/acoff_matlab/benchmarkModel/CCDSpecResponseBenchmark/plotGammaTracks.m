function plotGammaTracks(S)
%S is the track matrix. 
%Bijan Pourhamzeh - July 15th, 2008

%%%%%%%%%%%%New Matrix%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Column 1: Track ID
% Column 2: Step Number
% Column 3: Charge
% Column 4,5,6: Initial position (x,y,z)
% Column 7,8,9: Final position (x,y,z)
% Column 10: Track Length
% Column 11: Step Length
% Column 12: Final Energy (E)
% Column 13: dE

F=findstart(S); %Beginning of tracks
tracks=length(F);
hold on
%Plotting the tracks one by one (so that they are connected with a line)
for i=1:tracks-1
    plot3(S(F(i):F(i+1)-1,4),S(F(i):F(i+1)-1,5),S(F(i):F(i+1)-1,6))
end

%this is for the last track which isn't accounted for in the loop
plot3(S(F(length(F)):end,4),S(F(length(F)):end,5),S(F(length(F)):end,6))

%this plots a red dot at the end of each track
for j=1:tracks-1
    plot3(S(F(j+1)-1,4),S(F(j+1)-1,5),S(F(j+1)-1,6),'r.')
end

%this is for the last track
plot3(S(end,4),S(end,5),S(end,6),'r.')

end