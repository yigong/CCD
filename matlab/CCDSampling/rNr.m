function [x,y,phi]=rNr(xP,yP,xS,yS,pxsize)
%rotate and randomize starting position within a pixel of pitch, pxsize

randx=pxsize*(rand-0.5);
randy=pxsize*(rand-0.5);

phi=90*rand;
xPT=cosd(phi)*xP+sind(phi)*yP;
yP=-sind(phi)*xP+cosd(phi)*yP;
xP=xPT;
xST=cosd(phi)*xS+sind(phi)*yS;
yS=-sind(phi)*xS+cosd(phi)*yS;
xS=xST;

xP=xP+randx;
yP=yP+randy;
xS=xS+randx;
yS=yS+randy;


x=[xP; xS];
y=[yP; yS];