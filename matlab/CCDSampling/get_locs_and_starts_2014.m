function [no_primary_locs, no_secondary_locs, startpositions]=get_locs_and_starts_2014(in)
%dealing with omitting Edep=0 steps
%%How to check it worked: 
%>>[PRi,SCi,StartPos]=get_locs_and_starts_2014(**in**);
%>>**in**(StartPos(2)-2:StartPos(2)+6,2:end)


p_indices=find(in(:,11)==0);  %index of primaries = 0, " 2ndaries = 1;
dp=p_indices(2:end)-p_indices(1:end-1);  %length (in file lines) between primary indices
tmpp=find(dp~=1);
%tmp=find(in(:,5)==in(1,5));  %find where E=E_initialization... might have issues with a dE=0 in first step.
nhistories=length(tmpp)+1;

no_primary_locs=zeros(nhistories,1);
no_primary_locs(1)=tmpp(1);
for i=2:nhistories-1
  no_primary_locs(i)=tmpp(i)-tmpp(i-1);
end

%no_primary_locs=no_primary_locs-1;
no_primary_locs(end)=length(p_indices)-tmpp(end);


no_secondary_locs=zeros(nhistories,1);
no_secondary_locs(1)=p_indices(tmpp(1)+1)-no_primary_locs(1)-1;
for i=2:nhistories-1
    %no_secondary_locs(i)=p_indices(tmp(i))-no_primary_locs(i)-p_indices(tmp(i-1));
    no_secondary_locs(i)=p_indices(tmpp(i)+1)-no_primary_locs(i)-p_indices(tmpp(i-1)+1);
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


