function psft = MakePsfTable(tdt,D,xres,zres,savename)
%function psft = MakePsfTable(tdt,D,xres,zres,savename)
%
%Create a table of point-spread functions for charge carrier diffusion.
%Assumes 650 um thickness and 60 um width of psf kernel.
%
% ***Z coordinate reversed to measure from PIXEL PLANE not back plane
% (drift tables measured z from back plane.)
%
%tdt:   total drift time table in seconds (from TotalDriftTimes.mat)
%D:     diffusion coefficient (um^2/s)          (3.188e9)
%xres:  xy resolution of grid to use (um)       (0.5)
%zres:  z resolution of table, in um            (0.5)
%savename: filename to save table into.
%

z = 0:zres:650;
dt = interp1(tdt(:,1),tdt(:,2),(650-z'));   %z measured from pixel plane now

sigma = sqrt(2*D.*dt);  %um

gs = 60/xres + 1;    %gaussian size. should be an odd number. final units are xres.

psft = cell(size(sigma));
%create PSFs for each energy deposition point. [CPU intensive]
for i=1:length(sigma)
    psft{i} = fspecial('gaussian', [gs,gs],sigma(i,1)/xres);
end

%save
save(savename,'z','zres','psft')