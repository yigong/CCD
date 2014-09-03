function [no_primary_locs, no_secondary_locs, startpositions]=set_locs_and_starts(in,index)


p_indices=find(in(:,11)==index);  %index of primaries = 0, " 2ndaries = 1;
dp=in(p_indices(2:end),12)-in(p_indices(1:end-1),12);  %length (in file lines) between primary indices
tmp=(find(dp>1));  %finds where we skip between primaries
nhistories=length(tmp)+1;

no_primary_locs=zeros(nhistories,1);
no_primary_locs(1)=tmp(1);
for i=2:nhistories-1
  no_primary_locs(i)=tmp(i)-tmp(i-1);
end

%no_primary_locs=no_primary_locs-1;
no_primary_locs(end)=length(p_indices)-tmp(end);


no_secondary_locs=zeros(nhistories,1);
no_secondary_locs(1)=p_indices(tmp(1)+1)-no_primary_locs(1)-1;
for i=2:nhistories-1
    %no_secondary_locs(i)=p_indices(tmp(i))-no_primary_locs(i)-p_indices(tmp(i-1));
    no_secondary_locs(i)=p_indices(tmp(i)+1)-no_primary_locs(i)-p_indices(tmp(i-1)+1);
end
%no_secondary_locs(end-1)=p_indices(end)-no_primary_locs(end)-p_indices(tmp(end-1));
no_secondary_locs(end)=length(in(:,1))-p_indices(end);




startpositions=zeros(nhistories,1);
startpositions(1)=1;
for i=1:nhistories-1
    startpositions(i+1)=startpositions(i)+no_primary_locs(i)+no_secondary_locs(i);
end
%no_secondary_locs(end-1)=p_indices(end)-no_primary_locs(end-1)-p_indices(tmp(end));
%no_secondary_locs(end)=length(in(:,1))-no_primary_locs(end);


