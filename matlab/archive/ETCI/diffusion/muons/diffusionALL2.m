function [out]=diffusionALL(M,tdt,zp,theta,phi,plotflag)
%function [out]=diffusionALL(M,tdt,zp,theta,phi,plotflag)
%Applies diffusion015.m to all tracks in M. 
%Output is a cell array of SGrid for each track.

%Index each track
F=find(M(:,4)==0 & M(:,5)==0 & M(:,6)==0);

%Loop through all tracks except the last one
% progressbar2(0,0,'Diffusion');
%   ^-- progressbar doesn't work here, although it works anywhere else.
%       very mysterious. probably related to UFOs.
% progressfile;
out = cell(1,1000);
for i=1:length(F)-1
    out{i}=DiffuseTrack(M(F(i):F(i+1)-1,:),tdt,zp,theta,phi,plotflag);
%     progressbar2(i/length(F),0,'Diffusion');
%     progressfile(i/length(F));
end

%Last track
out{length(F)}=DiffuseTrack(M(F(length(F)-1):F(end),:),tdt,zp,theta,phi,plotflag);
% progressbar2(1,0,'Diffusion');
% progressfile(1);