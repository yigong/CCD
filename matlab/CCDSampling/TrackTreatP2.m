function [alphaO,phiO,dEO,hw]=TrackTreatP2(in)
%function [a,p,dE,hwt]=TrackTreatP2(in)
%is in index of tracks to look at, in is name of matlab parsed geant sim. always do both BP&PxP, backplane=(:,1) px plane=(:,2)


[no_primary_locs, no_secondary_locs, startpositions]=get_locs_and_starts_2014(in);
NL=length(no_primary_locs);
matlabpool;
%matlabpool size
NP=6;
NLpP=ceil(NL/NP);

idx=zeros(NP,2);
idx(1,1)=1; %start positions in the 'in' array for the next history
['segmenting input file: ',num2str(NP),' processors assinged']
hnos=zeros(NP,2);
hnos(1,1)=1;  %history #
for i=1:NP-1

    idx(i,2)=startpositions(NLpP*i)-1;
    in2{i}=in(idx(i,1):idx(i,2),:);
    idx(i+1,1)=startpositions(NLpP*i);
    hnos(i,2)=NLpP*i;
    hnos(i,1)=NLpP*(i-1)+1;
end
idx(NP,2)=length(in);
in2{NP}=in(idx(NP,1):idx(NP,2),:);
hnos(end,1)=NLpP*(NP-1)+1;
hnos(end,2)=NL;
clear in;

alphaO=zeros(NL,2);
phiO=zeros(NL,1);
dEO=phiO;
hw=alphaO;


parfor i=1:NP
    %f(i)=parfeval(runF,2,in2{i},i);    
    [a{i},p{i},dE{i},hwt{i}]=runF([hnos(i,1):hnos(i,2)],in2{i},i);
    pause(2)
end
for i=1:NP
    [alphaO(hnos(i,1):hnos(i,2),:)]=a{i};
    [phiO(hnos(i,1):hnos(i,2))]=p{i};
    [dEO(hnos(i,1):hnos(i,2))]=dE{i};
    [hw(hnos(i,1):hnos(i,2),:)]=hwt{i};
end
% for i=1:NP
%     [completedIdx,aR,pR,eR,hwR]=fetchNext(f);
%     alphaO(idx(completedIdx,1):idx(completedIdx,2),:)=aR;
%     phiO(idx(completedIdx,1):idx(completedIdx,2))=pR;
%     dEO(idx(completedIdx,1):idx(completedIdx,2))=eR;
%     hw(idx(completedIdx,1):idx(completedIdx,2),:)=hwR;
% end
matlabpool close;
end

function [alphaO,phiO,dEO,hw]=runF(is,in,PNo)

persistent psft 
load diffusioncode/psft.mat

rng(sum(clock)*1e3)

N=length(is);
alphaO=zeros(N,2);
phiO=zeros(N,1);
dEO=phiO;
hw=alphaO;
[no_primary_locs, no_secondary_locs, startpositions]=get_locs_and_starts_2014(in);
N=length(no_primary_locs);
for i=1:N
    [PNo i] 
    [P,S]=getHistory(in,i);
    [x,y,phi]=rNr(P(:,2), P(:,3), S(:,2),S(:,3),10.5);
    z=[P(:,4); S(:,4)];
    E=[P(:,6); S(:,6)];
    cd diffusioncode
    D=DiffuseTrack7(psft,'manual',[x y z E],'depthcoordinates',[.325 -.325]);  %%bk plane irradiation
    Dp=DiffuseTrack7(psft,'manual',[x y z E],'depthcoordinates',[-.325 .325]);  %%px plane irradiation       
    cd ..
    test=HybridTrack(D.img);
    testp=HybridTrack(Dp.img);
    if ~isfield(testp,'EdgeSegments');
        if ~isfield(test,'EdgeSegments');
            'no EdgeSegments'
            continue;
        end
        [test,hadwrong]=checkTrack(test,D);  %checks if radial position of 2nd found end is closer to zero than that of the first, if so, swaps them and retruns 1 for hadwrong. also returns -# if #ends~=2
        'no PxP EdgeSegments'
        if isfield(test,'alpha')
        alphaO(i,1)=mod(90-test.alpha,360);
        else
        alphaO(i,1)=NaN;
        end
        alphaO(i,2)=NaN;
        phiO(i)=phi;
        dEO(i)=sum(E);
        hw(i,1)=hadwrong;
        hw(i,2)=NaN;        
    elseif ~isfield(test,'EdgeSegments');
        [test,hadwrong]=checkTrack(testp,Dp);    %checks if radial position of 2nd found end is closer to zero than that of the first, if so, swaps them and retruns 1 for hadwrong. also returns -# if #ends~=2
        'no BP EdgeSegments'
        alphaO(i,1)=NaN;
        if isfield(testp,'alpha')
        alphaO(i,2)=mod(90-testp.alpha,360);
        else
        alphaO(i,2)=NaN;
        end
        phiO(i)=phi;
        dEO(i)=sum(E);
        hw(i,1)=NaN;
        hw(i,2)=hadwrong2;
    else        
        [test,hadwrong]=checkTrack(test,D);  
        [testp,hadwrong2]=checkTrack(testp,Dp);
        if isfield(test,'alpha')
        alphaO(i,1)=mod(90-test.alpha,360);
        else
        alphaO(i,1)=NaN;
        end
        if isfield(testp,'alpha')
        alphaO(i,2)=mod(90-testp.alpha,360);
        else
        alphaO(i,2)=NaN;
        end
        phiO(i)=phi;
        dEO(i)=sum(E);
        hw(i,1)=hadwrong;
        hw(i,2)=hadwrong2;
    end
end
end

function [test,hadwrong]=checkTrack(test,D)
%checks if radial position of 2nd found end is closer to zero than that of the first,
%If so, swaps them and retruns 1 for hadwrong. also returns -# if #ends~=2

ci=test.EdgeSegments.chosenIndex-1;  %should be 1 or two
if test.ends==2
    hadwrong=0;
    if sqrt((test.EdgeSegments.coordinatesPix(ci+1,1)+D.offsets(1)).^2+(test.EdgeSegments.coordinatesPix(ci+1,2)+D.offsets(2))^2)>sqrt((test.EdgeSegments.coordinatesPix(~ci+1,1)+D.offsets(1)).^2+(test.EdgeSegments.coordinatesPix(~ci+1,2)+D.offsets(2))^2)
        test=HybridTrackNO(D.img);
        hadwrong=1;
    end
else
    hadwrong=-test.ends;
end
end
