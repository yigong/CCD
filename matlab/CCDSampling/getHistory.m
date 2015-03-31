function [outP,outS]=getHistory(in,histno)

[pr,sc]=get_locs_and_starts_2014(in);
indx=1;
for i=1:histno-1
    indx=pr(i)+sc(i)+indx;
end
outP=in(indx:indx+pr(histno)-1,:);
outS=in(indx+pr(histno):indx+pr(histno)+sc(histno)-1,:);
