function [out]=diffusionALL(M,zp,theta,phi)
%Applies diffusion015.m to all tracks in M. 
%Output is a cell array of SGrid for each track.

%Index each track
F=find(M(:,5)==0 & M(:,6)==0 & M(:,7)==0);

%Loop through all tracks except the last one
for i=1:length(F)-1
    out{i}=diffusion015(M(F(i):F(i+1)-1,:),zp,theta,phi);
end

%Last track
out{length(F)}=diffusion015(M(F(length(F)-1):F(end),:),zp,theta,phi);
end