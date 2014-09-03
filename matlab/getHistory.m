function [outP,outS]=getHistory(in,histno)

[pr,sc]=set_locs_and_starts(in,0);
indx=1;
for i=1:histno-1
    indx=pr(i)+sc(i)+indx;
end
outP=in(indx:indx+pr(histno)-1,:);
outS=in(indx+pr(histno):indx+pr(histno)+sc(histno)-1,:);
