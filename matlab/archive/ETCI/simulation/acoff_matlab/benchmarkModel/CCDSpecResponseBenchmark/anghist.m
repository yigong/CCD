function [phi,ind]=anghist(S)
%Plots the angular distribution of compton scattered electrons.

%Only take electrons
S=S(S(:,3)==-1,:);
%Only take electrons which originate from original gamma (TID 1)
S=S(S(:,1)==2,:);
%Only take electron tracks with length greater than zero
F=findlongtracks(S);
L=length(F);
%Origin
O=[-1000.33*ones(L,1), zeros(L,1), zeros(L,1)];
%First and second step of electron track
pos0=S(F,4:6);
%pos1=S(F+1,4:6); % Old method of calculating e-trajectory
pos1=S(F,7:9);      % New way of calculating e-trajectory

rvec=pos0-O; %Vector from origin to first step
evec=pos1-pos0; %Vector from first to second step (initial direction)

%Magnitude of vectors
rvecnorm=sqrt(sum(rvec.*rvec,2));
evecnorm=sqrt(sum(evec.*evec,2));

%Dot product of vectors
dot=sum(rvec.*evec,2);

phi=acos(dot./(rvecnorm.*evecnorm))*180/pi;

BigTracks=F(phi>=90); %Tracks with scatter angle greater than 90 deg
ind=[];

%Find index of BigTracks
for i=1:length(BigTracks)
    ind=[ind; find(F==BigTracks(i))];
end

% for i=1:length(ind)-1
%     plot3(S(F(ind):F(ind+1)-1,4),S(F(ind):F(ind+1)-1,5),S(F(ind):F(ind+1)-1,6),'b.')
%     pause
% end

[n,x]=hist(phi,1000);
bar(x,n/L)
xlabel('Angle [deg]')
ylabel('Number of Electrons')
title('Angle Between Incident Gamma and Scattered Electron')

end


