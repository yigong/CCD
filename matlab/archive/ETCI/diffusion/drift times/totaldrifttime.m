function tdt = totaldrifttime(dt)
%function tdt = totaldrifttime(dt)
%
% Integrate differential drift time  from each depth,
%   to get total drift time

% load('DriftTime.mat')

%--CONTENTS OF dt--
%Col 1: depth z                    (0 = backside, 650 = pixels)    um
%Col 2: electric field E(z)        E(z) = [see below]               V/cm
%Col 3: drift velocity Vd          Vd(E) = interp1(Vd140)          cm/s
%Col 4: drift velocity Vd                                          um/s
%Col 5: differential drift time Td  = 1/Vd * dz                  s

tdt = zeros(size(dt,1),2);

tdt(1,1) = 650;
tdt(1,2) = dt(size(dt,1),5);

for i = 2:size(dt,1)
    tdt(i,1) = 650 - dt(i,1);
    tdt(i,2) = tdt(i-1,2) + dt(size(dt,1)+1-i,5);
end
% 
% [r,c] = find(dt(:,1)==z);  %assume z is part of the table
% 
% 
% 
% tdt = sum(dt(r:size(dt,1),4));