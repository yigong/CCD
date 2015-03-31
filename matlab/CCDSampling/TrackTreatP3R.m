function [alphaO,a,dEO,hw,betaO,ridgeWidth]=TrackTreatP3R(in,a,NTT,n)

%rotate input by a, divided into NTT serial calls n is index of serial call. always do both BP&PxP, backplane=(:,1) px plane=(:,2)


[no_primary_locs, no_secondary_locs, startpositions]=get_locs_and_starts_2014(in);
%%%%% NTT total tracktreatCalls, this is n of NTT
NL=length(no_primary_locs);
nh=ceil(NL/NTT);
if n<NTT
in=in(startpositions((n-1)*nh+1):startpositions(nh*n+1)-1,:);  %1: last line in 'in' before (nh*n+1) history
else
in=in(startpositions((n-1)*nh+1):end,:);
end

[no_primary_locs, no_secondary_locs, startpositions]=get_locs_and_starts_2014(in);
NL=length(no_primary_locs);
matlabpool;
%matlabpool size
NP=3;
NLpP=ceil(NL/NP);

idx=zeros(NP,2);
idx(1,1)=1; %start positions in the 'in' array for the next history
['segmenting input file: ',num2str(NP),' processors assinged']
hnos=zeros(NP,2);
hnos(1,1)=1;  %history #
for i=1:NP-1

    idx(i,2)=startpositions(NLpP*i+1)-1;
    in2{i}=in(idx(i,1):idx(i,2),:);
    idx(i+1,1)=startpositions(NLpP*i+1);
    hnos(i,2)=NLpP*i;
    hnos(i,1)=NLpP*(i-1)+1;
end
idx(NP,2)=length(in);
in2{NP}=in(idx(NP,1):idx(NP,2),:);
hnos(end,1)=NLpP*(NP-1)+1;
hnos(end,2)=NL; %NLpP*NP;
clear in;

alphaO=zeros(NL,2);
betaO=zeros(NL,2);
ridgeWidth=cell(NL,2);
dEO=zeros(NL,2);
hw=alphaO;


parfor i=1:NP
    %f(i)=parfeval(runF,2,in2{i},i);    
    [ao{i},bo{i},rwo{i},dE{i},hwt{i}]=runF([hnos(i,1):hnos(i,2)],in2{i},i,a);
    %[ao{i},p{i},dE{i},hwt{i}]=runF([hnos(i,1):hnos(i,2)],in2{i},i,a);
    pause(1); %pause for 2s to ensure random#sead in runF is randomized
end
for i=1:NP
    length(ao{i})
    hnos(i,2)-hnos(i,1)+1
    [alphaO(hnos(i,1):hnos(i,2),:)]=ao{i};
    [betaO(hnos(i,1):hnos(i,2),:)]=bo{i};
    [ridgeWidth(hnos(i,1):hnos(i,2),:)]=rwo{i};
    %[phiO(hnos(i,1):hnos(i,2))]=p{i};
    [dEO(hnos(i,1):hnos(i,2),:)]=dE{i};
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

function [alphaO,betaO,ridgeWidth,dEO,hw]=runF(is,in,PNo,a)
%function [alphaO,phiO,dEO,hw]=runF(is,in,PNo,a)

persistent psft 
load diffusioncode/psft.mat

rng(sum(clock*1e3));

N=length(is);
alphaO=zeros(N,2);
betaO=zeros(N,2);
dEO=zeros(N,2);
hw=alphaO;
ridgeWidth = cell(N,2);


[no_primary_locs, no_secondary_locs, startpositions]=get_locs_and_starts_2014(in);
N=length(no_primary_locs);
for i=1:N
    [PNo i] 
    [P,S]=getHistory(in,i);
    [x,y]=rnrF(P(:,2), P(:,3), S(:,2),S(:,3),10.5,a);  %phi=a
    z=[P(:,4); S(:,4)];
    E=[P(:,6); S(:,6)];
    cd diffusioncode
    D=DiffuseTrack7(psft,'manual',[x y z E],'depthcoordinates',[.325 -.325]);  %%bk plane irradiation
    Dp=DiffuseTrack7(psft,'manual',[x y z E],'depthcoordinates',[-.325 .325]);  %%px plane irradiation       
    cd ..
    test=HybridTrack(D.img);
    testp=HybridTrack(Dp.img);
    if ~isfield(test,'err') %bk plane irradiation
        [test,hadwrong]=checkTrack(test,D);  %checks if radial position of 2nd found end is closer to zero than that of the first, if so, swaps them and retruns 1 for hadwrong. also returns -# if #ends~=2
        %'no PxP EdgeSegments'
        if isfield(test,'alpha') %alpha
        alphaO(i,1)=mod(90-test.alpha,360);
        else
        alphaO(i,1)=NaN;
        end
        if isfield(test, 'beta') %beta
        betaO(i,1) = test.beta;
        else
        betaO(i,1) = NaN;
        end
        if isfield(test, 'w') %ridgeWidth
        ridgeWidth{i,1} = test.w;
        else
        ridgeWidth{i,1} = NaN;
        end
        dEO(i,1)=test.Etot;
        hw(i,1)=hadwrong;
    else
        alphaO(i,1)=NaN;
        betaO(i,1) = NaN;
        ridgeWidth{i,1} = NaN;
        dEO(i,1)=test.Etot;
        hw(i,1)=hadwrong;
    end
    
    if ~isfield(testp,'err') % px plane irradiation
        [testp,hadwrong]=checkTrack(testp,Dp);    %checks if radial position of 2nd found end is closer to zero than that of the first, if so, swaps them and retruns 1 for hadwrong. also returns -# if #ends~=2
        if isfield(testp,'alpha') %alpha
        alphaO(i,2)=mod(90-testp.alpha,360);
        else
        alphaO(i,2)=NaN;
        end        
        if isfield(testp, 'beta') %beta
        betaO(i,2) = testp.beta;
        else
        betaO(i,2) = NaN;
        end        
        if isfield(testp, 'w') %ridgeWidth
        ridgeWidth{i,2} = testp.w;
        else
        ridgeWidth{i,2} = NaN;
        end
        dEO(i,2)=testp.Etot;
        hw(i,2)=hadwrong;
    else
        alphaO(i,2)=NaN;
        betaO(i,2) = NaN;
        ridgeWidth{i,2} = NaN;
        dEO(i,2)=test.Etot;
        hw(i,2)=hadwrong;

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
