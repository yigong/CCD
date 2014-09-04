function [Etot,Edep,Eesc,multiplicity] = GetEnergiesFromGeant(M)
%function [Etot,Edep,Eesc,multiplicity] = GetEnergiesFromGeant(M)
%
% Based on DiffuseTrack3, 4/22/11.
%
% New in DiffuseTrack3:
%  geant4 must write any particle step that produces a secondary.
%    This means that secondary electrons need not be matched to their parent,
%    because there is no ambiguity in deposited energy.
%  multiple electrons (from a single photon parent) are all tracked, and output
%    as cells.
%
%input variables
%
%M: electron data matrix from geant4. Columns: (units)
%       1: Track ID
%       2: Parent ID
%       3: Step num
%       4: Charge                           (e-)
%       5,6,7: Initial position (x,y,z)     (mm)
%       8,9,10: Final position (x,y,z)      (mm)
%       11: Accumulated track length        (mm)
%       12: Step length                     (mm)
%       13: Final energy (after step)       (eV)
%       14: Initial energy (before step)    (eV)
%       15: Energy deposited (dE)           (eV)
%psft: table of point spread functions, from MakePsfTable(tdt,D,0.5,zres,savename)
%      cell array of 121x121 psf's
%th: thickness of CCD           (um)
%
% Coordinate system: x_geant is depth dimension (pixel face = 0, back face = (th/1000))
%
%[blsig: standard deviation of black level (background noise), in keV.]
% [blsig set to 0.01874 keV, as measured on 11/2007 data.]
%[Fano factor assumed to be 0.14]
%
%plotflag: 1 to plot CCD image
%          0 for no plot
%
%output variables
%
%D is a cell array, with one cell for each primary electron in this event sequence.
%
%D{i}.img:     data matrix of ccd plane pixels.    (keV)
%D{i}.cheat:
%     cheat.Etot:    actual electron energy        (keV)
%     cheat.Edep:    energy deposited in detector  (keV)
%     cheat.a:       alpha                         (degrees)
%     cheat.b:       beta                          (degrees)
%     cheat.z:       zp (starting z)               (um)
%     cheat.contained: logical for all energy deposited in active volume
%                  (due to rounding errors, Edep ~= Etot)
%     cheat.x:     x coordinate of start of trajectory (pixels)
%     cheat.y:     y coordinate of start of trajectory (pixels)
%D{i}.x2:      list of 3D coordinates that have been diffused. z=depth.
%D{i}.dE:      list of energies diffused. corresponds to D.x2.

% pixsize = 10.5;     %um
% buffer = 5.5*pixsize;    %um            %NOT ANY LONGER: point (0,0) WAS in the center of a pixel

% zres = 0.5; %um resolution of PSFT
% zres = th/(length(psft)-1);
% int=0.5; %um    %int must match the xres of PSFT
%This is the resolution of the fine grid. Energy deposition coordinates
%will be rounded off to a grid position.

%shortcuts to column numbers
indTrackID = 1;
indParentID = 2;
indStepNum = 3;
indCharge = 4;
indInitPos = 5:7;
indFinalPos = 8:10;
indTrackLen = 11;
indStepLen = 12;
indFinalE = 13;
indInitE = 14;
inddE = 15;

%~~~~
%start rewritten GetEnergiesFromGeant code

N = max(M(:,indTrackID));

f1 = nan(1,N);
f2 = nan(1,N);
E1 = nan(1,N);
E2 = nan(1,N);
p = nan(1,N);

for i=1:N
    fcur = find(M(:,indTrackID)==i,1,'first');
    if ~isempty(fcur)
        f1(i) = fcur;
        f2(i) = find(M(:,indTrackID)==i,1,'last');
        E1(i) = M(f1(i),indInitE);  %particle start energy
        E2(i) = M(f2(i),indFinalE); %particle escape energy
        p(i) = M(f1(i),indParentID);
        q(i) = M(f1(i),indCharge);
    else
        f1(i) = nan;
        f2(i) = nan;
        E1(i) = nan;
        E2(i) = nan;
        p(i) = nan;
        q(i) = nan;
    end
end

primaryelectrons = find(p==1 & q==-1);
multiplicity = length(primaryelectrons);

pg(1:N) = nan;
pg(p==1) = find(p==1);
while true
    pgtmp = pg;
    for i=1:length(pg)
        if isnan(pg(i)) && ~isnan(p(i)) && p(i)>0
            pg(i) = pg(p(i));
        end
    end
    if all( (pgtmp==pg) | (isnan(pgtmp) & isnan(pg)))
        break
    end
end

Etot = nan(1,multiplicity);
Edep = nan(1,multiplicity);
Eesc = nan(1,multiplicity);

for i=1:length(primaryelectrons)
    Etot(i) = E1(primaryelectrons(i));
    Eesc(i) = sum(E2(pg==primaryelectrons(i)));
end
Edep = Etot - Eesc;

Edep = Edep / 1e3;  %eV to keV
Etot = Etot / 1e3;
Eesc = Eesc / 1e3;

return


%end rewritten GetEnergiesFromGeant code
%~~~~~~~

%% Particle info.
%just some extra shortcuts to use later.

ind = 1;
% f1 = zeros(1,N);
% % f2 = zeros(1,N);
% q = zeros(1,N);
% Ep = zeros(1,N);
% p = zeros(1,N);

asdf=0; %dummy variable for Matlab 2008b

% contained = true(1,N);   %until proven otherwise... see line 344ff, this is "new"
% Eescaped = nan(1,N);
for i=1:max(M(:,indTrackID))
    fcur = find(M(:,indTrackID)==i,1,'first');  %first record of this particle
    if ~isempty(fcur)
        id(ind) = i;
        id2ind(i) = ind;
        f1(ind) = fcur;               %first record of this particle
    %     f2(i) = find(M(:,1)==i,1,'last');   %last record of this particle
        q(ind) = M(fcur,indCharge);          %charge
        Ep(ind) = M(fcur,indInitE);        %initial energy
        if ind>1
            if ~isempty(find(id==M(fcur,indParentID)))
                p(ind) = find(id==M(fcur,indParentID));          %parent ID index
            else
                id(ind) = [];
                f1(ind) = [];
                q(ind) = [];
                Ep(ind) = [];
                continue
            end
        else
            p(ind) = 0;
        end
        %check for escape
        Eescaped(ind) = M(find(M(:,indTrackID)==i,1,'last'),indFinalE);
        if Eescaped(ind)>0
            contained(ind) = false;
        else
            contained(ind) = true;
        end    
        ind = ind+1;
    end
end
N = ind-1;    %particles tracked

%% Count primary electrons.

%particle #1 should be the parent photon
if q(1) > 0 || q(1) < 0     %allow for the possibility of photon not being recorded,
                            %  in which case q(1) = nan
    error('Particle #1 is not a photon..')
end

%particlestatus: added 4/21/11.
% 1 = primary photon (defined as particle ID=1)
% 2 = primary electron (parentID = 1 and charge = -1)
% 3 = secondary, etc. electron (charge = -1 and charge(parentID) = -1)
% 4 = secondary (etc) photon (bremsstrahlung) (charge = 0 and charge(parentID) = -1)
% 5 = brems-induced electron (charge = -1 and charge(parentID) = 0 and parentID >1)
% 6 = brems-induced secondary electron (charge = -1 and charge(parentID) = -1 and particlestatus(parentID) = 5)

if id(1)~=1
    error('Particle #1 is not in the first index...')
    %this is assumed below...
end
particlestatus(1:N) = nan;
particlestatus(id==1) = 1;                              %primary photon
particlestatus(p==1 & q==-1) = 2;                %primary electron
particlestatus(q==-1 & [false,q(p(2:end))==-1]) = 3;    %secondary electron
particlestatus(q==0 & [false,q(p(2:end))==-1]) = 4;     %secondary photon
particlestatus(q==-1 & [false,q(p(2:end))==0] & p>1) = 5; %brems-induced electron

% %secondary brems-induced electrons
% while true
%     pstmp = particlestatus;
%     particlestatus(particlestatus==3 & (particlestatus(p)==5 | particlestatus(p)==6)) = 6;
%     if all(pstmp==particlestatus)
%         break
%     end
% end

multiplicity = sum(particlestatus==2);  %number of primary electrons
bremsmultiplicity = sum(particlestatus==5); %number of brems-induced electrons

%separate all particles into groups of contiguous deposition.
particlegroup(1:N) = nan;
particlegroup(particlestatus==2) = find(particlestatus==2);     %primary electrons keep their index number
particlegroup(particlestatus==4) = find(particlestatus==4);     %brems photons keep their index number
%particles inherent particlegroup from their parent.
while true
    pgtmp = particlegroup;
    particlegroup(isnan(particlegroup) & ~isnan([nan,particlegroup(p(2:end))])) = particlegroup(p(isnan(particlegroup) & ~isnan([nan,particlegroup(p(2:end))])));
    if all(pgtmp==particlegroup | (isnan(pgtmp) & isnan(particlegroup)))
        break
    end
end
%number particlegroups from 1 to #
%also: identify as brems or primary
%also: translate particle containment, Eescaped to particlegroup containment and Eescaped.
[upg,pgind,asdf] = unique(particlegroup,'first');
bcontained(1:multiplicity) = true;  %all brems from this particle is contained
pcontained(1:multiplicity) = true;  %all primary electron deposition from this particle is contained
pEescaped(1:multiplicity) = 0;      %total escaped energy from this primary electron and daughters
for i=1:length(upg)
    if isnan(upg(i))
        continue    %skip primary photon
    end
    bremsflag(i) = particlestatus(pgind(i))==5;    %true for brems groups, false for primary electron groups
    if bremsflag(i) && ~all(contained(particlegroup==upg(i)))
        %index of associated primary electron
        pind = find(upg==particlegroup(p(upg(i))));
        bcontained(pind - sum(bremsflag(1:pind-1))) = false;
        pEescaped(pind - sum(bremsflag(1:pind-1))) = pEescaped(pind - sum(bremsflag(1:pind-1))) + sum(Eescaped(particlegroup==upg(i)));
    else
        pcontained(i - sum(bremsflag(1:i-1))) = all(contained(particlegroup==upg(i)));
        pEescaped(i - sum(bremsflag(1:i-1))) = pEescaped(i - sum(bremsflag(1:i-1))) + sum(Eescaped(particlegroup==upg(i)));
    end
    
    particlegroup(particlegroup==upg(i)) = i;
end

%this will be needed in marking pileup stuff later.

% %above section replaced:
%{
% %% Remove multiple Compton interactions.
% % Only look at first Compton electron.
% fCompton = find(p==1);  %expect #1 to be the photon...
% if length(fCompton)>1
%     [asdf,mC] = min(M(f1,1));
%     fCompton(mC) = [];  %remove the first Compton electron from this list.
%     
%     %remove the others from the matrix.
%     i=1;    %index for fCompton list
%     while i <= length(fCompton)
%         %flag daughters of this particle for removal
%         fadd = find(p==fCompton(i));
%         fCompton = [fCompton, fadd];
%         %remove this particle
%         M(M(:,1)==fCompton(i),:) = [];
%         i = i+1;
%     end
%     
%     %redefine some variables
%     N = max(M(:,1));
%     f1 = NaN(1,N);
%     for i=1:N
%         if ~isempty(find(M(:,1)==i,1,'first'))
%             f1(i) = find(M(:,1)==i,1,'first');
%         end
%     end
% end
%}
%% Step info.

%need to record dE (col 15) for every step written,
%  and delta-E (col14 - col13) in between consecutive steps of one particle.

%define variables: 
%   x (position of energy deposition, (mm,mm,mm))
%   dE (amount of energy deposited, eV)
%   ss (step size.. do I care?)

%initalize vectors with enough space for everything (two entries for every row of M)
% intra-step (within step i): index 2*i
% inter-step (between step i and i+1): index 2*i+1
x = nan(size(M,1)*2,3);
dE = nan(size(M,1)*2,1);
pg = nan(size(M,1)*2,1);
% ss = nan(size(M,1)*2,1);

%first: photons "deposit" the binding energy of the electron they impact.
%       this energy will result in an X-ray of <= 1.8 keV (in Si) which will deposit
%       essentially at the position of the interaction. (very short attenuation length)
fcur = find(M(:,indCharge)==0);    %assume no other particles tracked that have charge 0.
x(2*fcur,:) = M(fcur,indFinalPos);   %assign deposited energy to end of photon step, where the interaction is
dE(2*fcur,1) = M(fcur,inddE);
for i=1:length(fcur)
    if ~isnan(particlegroup(M(fcur(i),indTrackID)))    %brems photons belong to their own particle group.
        pg(2*fcur(i),1) = particlegroup(M(fcur(i),indTrackID));
    else    %primary photon does not belong to any particle group. bummer.
        matchf = find((x(2*fcur(i),1)==M(:,indInitPos(1))) & ...
            (x(2*fcur(i),2)==M(:,indInitPos(2))) & ...
            (x(2*fcur(i),3)==M(:,indInitPos(3)))); %i think the final photon position is exactly the initial electron position?
        if isempty(matchf)  %otherwise, get closest
            Mtmp = M(~isnan(particlegroup(id2ind(M(:,indTrackID)))),indInitPos);
            [~,indtmp] = min(norm2(repmat(x(2*fcur(i),1:3),size(Mtmp,1),1)-Mtmp));
            indtmp = find(~isnan(particlegroup(id2ind(M(:,indTrackID)))),indtmp,'first');
            matchf = indtmp(end);
        end
        pg(2*fcur(i),1) = particlegroup(id2ind(M(matchf,indTrackID)));
    end
end

%second: intra-step electron deposition.
fcur = find(M(:,indCharge)==-1);  %assume no other particles tracked that have charge -1. (If there are, they deposit energy too...)
x(2*fcur,:) = (M(fcur,indFinalPos) + M(fcur,indInitPos))./2; %average initial, final positions
dE(2*fcur,1) = M(fcur,inddE);
pg(2*fcur,1) = particlegroup(M(fcur,indTrackID));

%second: inter-step electron deposition. put this in the "second" position.
lgSameParticle = [(diff(M(:,indTrackID))==0); 0];   %lgSameParticle(i) = true if rows i, i+1 are from the same particle
fcur = find(lgSameParticle & M(:,indCharge)==-1);
x(2*fcur+1,:) = (M(fcur,indFinalPos) + M(fcur+1,indInitPos))./2; %average final position i and initial position i+1
dE(2*fcur+1,1) = M(fcur,indFinalE) - M(fcur+1,indInitE);
pg(2*fcur+1,1) = particlegroup(M(fcur,indTrackID));

%previous code:
%{
% %Inter-record steps.
% x = [0,0,0; (M(1:end-1,8:10) + M(2:end,8:10))./2];   %average of final positions from one record to next.
% ss = [0;norm2(diff(M(:,8:10)))];                %step sizes from one record to the next.
% dE = [0;-diff(M(:,13))];        %difference in final energies of each record.
% %A given entry contains the difference from the previous line, to this line.
% 
% %Intra-record steps (at beginning of each particle track)
% x(f1(~isnan(f1)),:) = (M(f1(~isnan(f1)),5:7) + M(f1(~isnan(f1)),8:10))./2;
% ss(f1(~isnan(f1))) = norm2(M(f1(~isnan(f1)),5:7) - M(f1(~isnan(f1)),8:10));
% dE(f1(~isnan(f1))) = M(f1(~isnan(f1)),14) - M(f1(~isnan(f1)),13);
% %note: I have used [14]-[13] instead of [15] because I am expecting to see the energy of the secondary in its parent.
% 
% dE(M(:,4)==0) = M(M(:,4)==0,15);    %photons: dE represents binding energy of Compton-/photo-electron, 
% x(M(:,4)==0,:) = M(M(:,4)==0,8:10);   %  which will be deposited by x-ray very near the atom.
%}
%% Deal with secondaries.
%must subtract energy carried by secondary particle, from energy lost by parent.
%xxxxx NO LONGER NEEDED!

%{
%Do this carefully to avoid misidentifying parent steps. 
%   There may be two or more steps that are close enough and have enough dE.

fparent = cell(1,N);
unfinished = false(size(fparent)); %will be set to true as we go

maxdistance = 1.2e-3;   %mm

Eapp = false(1,N);

%first loop: only apply match if there is a single solution.
for i=2:N           %first particle has no parent.
    if q(p(i))~=-1 ...      %if parent is not an electron (no need to subtract dE)
            || isnan(f1(i))  % or this particle was removed as a secondary Compton (dE does not exist anyway)
        %ignore.
        fparent{i} = [];
        continue
    end
    f = find(M(:,1)==p(i) ...    %ID matches parent ID given
        & dE >= Ep(i) ...         %energy loss is large enough
        & norm2(x - repmat(M(f1(i),5:7),length(dE),1)) < max(maxdistance*ones(length(dE),1),ss));  %less than 1.2 um or step size
    if isempty(f) && (M(f1(i),5) < 1e-3 || M(f1(i),5) > (th-1)*1e-3)
        %escape error: see file index 999 #1591
        %the parent particle energy deposition has not been tracked anyway.
        %so just let the secondary deposit its energy, no worries.
        fparent{i} = [];
    elseif isempty(f)
        disp(['Parent step not found for particle #',num2str(i),'!'])
        fparent{i} = [];
        %*************** error
    elseif length(f)>1
        fparent{i} = f;
        unfinished(i) = true;
        continue
    elseif length(f)==1
        %subtract energy of secondary particle now.
        dE(f) = dE(f) - Ep(i);
        Eapp(i) = true;
        fparent{i} = [];
    else
        %????
        disp('wth??')
    end
%     [m,ind] = min(norm2(x(fparent,:) - repmat(M(f1(i),5:7),length(fparent),1)));
%     
%     if m > max(1e-3,ss(fparent(ind)))   %distance is greater than step size, and greater than 1 um. step size calculated from positions, not "step size" column. (see above)
%         disp(['Parent step too far from secondary! secondary ID ',num2str(i),', distance = ',num2str(m),' mm'])
% %         continue
%     end
    
end

%second loop: 
changedflag = true;
guessflag = false;

while any(unfinished)    %exit loop when finished
    if ~changedflag
        guessflag = true;
    end
    changedflag = false;
    
    %first, make a guess if that is needed.
    if guessflag
        %guess on the highest energy secondary.
        [asdf,ig] = max(Ep(unfinished));
        funf = find(unfinished);
        ig = funf(ig);
        %take minimum distance match.
        [asdf,ind] = min(norm2(x(fparent{ig},:) - repmat(M(f1(ig),5:7),length(fparent{ig}),1)));
        dE(fparent{ig}(ind)) = dE(fparent{ig}(ind)) - Ep(ig);
        fparent{ig} = [];
        unfinished(ig) = false;
        Eapp(ig) = true;
        %now, check again for any solutions that this caused.
    end
    
    for i=2:N
        if ~isempty(fparent{i})
            fparent{i}(dE(fparent{i}) < Ep(i)) = [];    %clear matches that are no longer valid due to dE.
            
            if length(fparent{i})==1
                %subtract energy of secondary particle now.
                dE(fparent{i}) = dE(fparent{i}) - Ep(i);
                fparent{i} = [];
                changedflag = true;
                unfinished(i) = false;
                Eapp(i) = true;
            elseif length(fparent{i}) > 1
                %check for unique matches.
                u = true(size(fparent{i})); %matches that are unique to this secondary
                for j=[2:(i-1),(i+1):N]
                    [m1,m2] = meshgrid(fparent{i},fparent{j}');
                    if any(m1(:)==m2(:))
                        [asdf,mf] = find(m1==m2);
                        u(mf) = false;
                    end
                end
                if sum(u)==1    %one unique match.
                    %subtract energy of secondary particle now.
                    dE(fparent{i}(u)) = dE(fparent{i}(u)) - Ep(i);
                    fparent{i} = [];
                    changedflag = true;
                    unfinished(i) = false;
                    Eapp(i) = true;
                elseif sum(u) > 1   %several unique matches.
                    %take minimum distance match.
                    [asdf,ind] = min(norm2(x(fparent{i}(u),:) - repmat(M(f1(i),5:7),sum(u),1)));
                    uf = find(u);
                    dE(fparent{i}(uf(ind))) = dE(fparent{i}(uf(ind))) - Ep(i);
                    fparent{i} = [];
                    changedflag = true;
                    unfinished(i) = false;
                    Eapp(i) = true;
                end
            else
                %no matches available! (this could happen after guessing, perhaps)
                error(['Cannot find match for particle #',num2str(i),'!'])
            end
        end
    end
end

%}

%% Unit conversion, coordinate transform.

x = x .*1e3;        %mm to um
dE = dE .*1e-3;     %eV to keV

%new: dimensions stay the same... 4/21/11

%old:
% %  x(:,1) is the depth dimension. ('x' => 'z')
% %  x(:,2) is the first pixel dimension. ('y' => 'x')
% %  x(:,3) is the second pixel dimension. ('z' => 'y')
% x2(:,3) = x(:,1);
% x2(:,1) = x(:,2);
% x2(:,2) = x(:,3);
% x = x2; clear x2

%flag escaped electrons

%4/22/11: geant4 now discarding electron steps outside detector volume, so we can't
%use that to see escapes anymore. instead, look for particle final energy > 0
%i put this above, in the 1:N loop.

%old:
escapes = x(:,3) < -th | x(:,3) > 0;     %4/21/11: flipping signs so that detector is from 0 to -th
% % escapes = x(:,3) > th | x(:,3) < 0;
% 
% %'escapes' is twice as long as M; so index M doubly.
% if any(escapes & M(ceil(0.5:0.5:size(M,1)),4)==-1)
%     contained = false;
% else
%     contained = true;
% end

x = x(~escapes & dE>0,:);     
%without dE>0, empty energy deposition (e.g. from movement of source photon) might cause the track image to be huge
%this also gets rid of NaN's from the initialization of x and dE (e.g. from the space for inter-record deposition, where not applicable)
dE = dE(~escapes & dE>0);

%% pre-segment tracks if necessary
% to avoid memory allocation of a huge image.
% 0.5 um fine grid and 10.5 um pixels: >400 grid points per pixel.
% If two Compton events happen to be located at opposite corners of the CCD, then we have
%  1454 x 726 x 441 x 4 bytes = 1.8GB. Probably okay on server, but let's try to avoid that.

% %matrix for pileup of particlegroups.
% overlap = false(max(particlegroup),max(particlegroup));

%flags for pileup of primary electrons or of brems (to go in cheat structure)
pileup(1:multiplicity) = 0; %start with no pileup
brems(1:multiplicity) = 0;  %start with no brems
pileupind = 1;  %index for pileup groups

%create image from each particlegroup separately.
for i=1:max(particlegroup)
    %image bounds - rounded to next fine grid mark..
    minx(i) = int * floor((min(x(pg==i,1)) - buffer) / int);
    miny(i) = int * floor((min(x(pg==i,2)) - buffer) / int);
    maxx(i) = int * ceil((max(x(pg==i,1)) + buffer) / int);
    maxy(i) = int * ceil((max(x(pg==i,2)) + buffer) / int);
    for j=1:(i-1)
        if minx(j)<maxx(i) && minx(i)<maxx(j) && miny(j)<maxy(i) && miny(i)<maxy(j)
            %set pileup, brems markings
            if ~bremsflag(find(particlegroup==i,1)) && ...
                    ~bremsflag(find(particlegroup==j,1)) %two primary electrons
                %mark as a pileup group
                if pileup(j - sum(bremsflag(1:j-1))) > 0    %group j already marked as a pileup group
                    pileup(i - sum(bremsflag(1:i-1))) = pileup(j - sum(bremsflag(1:i-1))); %same group
                else
                    pileup(i - sum(bremsflag(1:i-1))) = pileupind;
                    pileup(j - sum(bremsflag(1:i-1))) = pileupind;
                    pileupind = pileupind + 1;
                end
            elseif ~bremsflag(find(particlegroup==i,1)) && ...
                    bremsflag(find(particlegroup==j,1)) %i primary and j brems
                %count another brems on this primary track
                brems(i - sum(bremsflag(1:i-1))) = brems(i - sum(bremsflag(1:i-1))) + 1;
            elseif bremsflag(find(particlegroup==i,1)) && ...
                    ~bremsflag(find(particlegroup==j,1)) %i brems and j primary
                %count another brems on this primary track
                brems(j - sum(bremsflag(1:j-1))) = brems(j - sum(bremsflag(1:j-1))) + 1;
            else    %two brems
                %don't care.
            end
            
            %now combine particlegroups.
            pg(pg==i) = j;
        end
    end
end

upg = unique(pg);   %unique images to generate
clear minx maxx miny maxy

gs = size(psft{1},1);    %gaussian size of psf.
gs2 = (gs-1)/2;         %radial size of psf

F = 0.14;           %Fano factor to use. Silicon: 0.08 < F < 0.14 (Knoll)
blsig = 0.01874;    %std dev of background level, in keV. Nov 2007 data: 0.01874
epcc = 3.7e-3;      %energy per eh pair, in keV
    

for i=1:length(upg)
    
    lgcur = pg==upg(i); %points to work with in this image
    xcur = x(lgcur);
    dEcur = dE(lgcur);
    
    %% Remove stray points - no longer necessary.
    %  (e.g. distant brehmsstrahlung deposition, ???)

    %old code.
    %{
    loopcount = 0;

    while true
        %While stray points are still being removed...

        %image bounds
        minx=min(x(:,1));
        miny=min(x(:,2));
        maxx=max(x(:,1));
        maxy=max(x(:,2));

        x2(:,1) = x(:,1) + (buffer-minx);  %microns from image edge.
        x2(:,2) = x(:,2) + (buffer-miny);  %microns from image edge.
        x2(:,3) = x(:,3);   %depth

        %position of start of track
        x_0 = x2(1,1) / pixsize;
        y_0 = x2(1,2) / pixsize;

        %imsize in pixels
        imsize = [ceil((maxx + 2*buffer-minx)/pixsize), ...
            ceil((maxy + 2*buffer-miny)/pixsize)];

        %check for stray points.
        if max(imsize) > 100    %arbitrary threshold..
            [nproj,xproj] = hist(x2(:,1),buffer:50:imsize(1)*pixsize-buffer);   %projection along x-axis
            if any(nproj==0)
                %gap in histogram: separation of one part of "track" from the rest
                %identify "main" part of track
                n1 = sum(nproj(1:find(nproj==0,1)));
                n2 = sum(nproj(find(nproj==0,1,'last'):end));
                if n1 > n2
                    %remove n2 end of image.
                    xt = xproj(find(nproj==0,1,'last'));
                    %flag for removal
                    rmlogic = x2(:,1) > xt;
                else
                    %remove n1 end of image.
                    xt = xproj(find(nproj==0,1));
                    %flag for removal
                    rmlogic = x2(:,1) < xt;
                end

                %remove
                x(rmlogic,:) = [];
                dE(rmlogic,:) = [];
                clear x2
            else
                %check projection along y-axis
                [nproj,xproj] = hist(x2(:,2),buffer:50:imsize(2)*pixsize-buffer);
                if any(nproj==0)
                    n1 = sum(nproj(1:find(nproj==0,1)));
                    n2 = sum(nproj(find(nproj==0,1,'last'):end));
                    if n1 > n2
                        %remove n2 side
                        xt = xproj(find(nproj==0,1,'last'));
                        rmlogic = x2(:,2) > xt;
                    else
                        %remove n1 side
                        xt = xproj(find(nproj==0,1));
                        rmlogic = x2(:,2) < xt;
                    end

                    %remove
                    x(rmlogic,:) = [];
                    dE(rmlogic,:) = [];
                    clear x2
                else
                    %both projections look good: proceed
                    break
                end
            end
        else
            break
        end
        loopcount = loopcount + 1;
        if loopcount > 5
            asdf = 0;
        end
    end

    %}

    %% Set up image

    %image bounds - rounded to next pixel...
    minx = pixsize * floor((min(xcur(:,1)) - buffer) / pixsize);
    miny = pixsize * floor((min(xcur(:,2)) - buffer) / pixsize);
    maxx = pixsize * ceil((max(xcur(:,1)) + buffer) / pixsize);
    maxy = pixsize * ceil((max(xcur(:,2)) + buffer) / pixsize);
    
    %position adjustment to image coordinates. still in um.
    x2(:,1) = xcur(:,1) - minx;
    x2(:,2) = xcur(:,2) - miny;
    x2(:,3) = xcur(:,3);
    
    %image size in pixels
    imsize = [(maxx-minx)/pixsize, (maxy-miny)/pixsize];

    %initialize fine grid of track image.
    fgrid = zeros(imsize.*pixsize./int);
    
    %lookup psf for each energy deposition point.
    %sum into image of resolution 'int'.
    for j=1:size(x2,1)
        zi = round(abs(x2(j,3))/zres)+1;     %index of psft to use. assuming 0 is pixel face
        Espread = dEcur(j,1).*psft{zi};    %(already keV units)
        xi = round(x2(j,1)/int);        %round x,y to the nearest fine grid coordinate
        yi = round(x2(j,2)/int);
        fgrid(xi-gs2:xi+gs2,yi-gs2:yi+gs2) = fgrid(xi-gs2:xi+gs2,yi-gs2:yi+gs2) + Espread;
    end

    %create CCD image at pixel resolution.
    %start as cell array, we will pull it back out of the cell at the end of it is only a single image.
    img{i}=zeros(imsize);
    for s=1:imsize(1)
        for r=1:imsize(2)
    %         SGrid(s,r)=sum(sum(Grid(51+pxsize*(s-1)/int:50+pxsize*(s)/int,51+pysize*(r-1)/int:50+pysize*(r)/int)));  
            img{i}(s,r) = sum(sum(fgrid(1+(s-1)*pixsize/int:s*pixsize/int,1+(r-1)*pixsize/int:r*pixsize/int)));
        end
    end

    %Noise: dark current + Fano
    cc = img{i} ./ epcc; %average # charge carriers in each pixel
    u = sqrt(blsig.^2 + epcc^2 .* F .* cc);    %add stdev's in quadrature to get total sigma of noise

    img{i} = img{i} + u.*SampleGauss2(size(img));     %add gaussian noise of magnitude u over image

    step1 = f1(find(q(particlegroup==i)==-1,1));
    %Set up outputs
    cheat.Edep(i) = sum(dEcur);     %already in keV
    cheat.Etot(i) = M(step1,indInitE) .* 1e-3; %first E_initial, in keV
    %due to rounding errors, abs(Etot - Edep) = ~4e-7 *Etot, for contained tracks

    cheat.x(i) = x2(1,1) / pixsize;  %in pixel coordinates of image
    cheat.y(i) = x2(1,2) / pixsize;  %in pixel coordinates of image
    cheat.z(i) = x2(1,3);  %already in um

    %initial direction
    dx = M(step1,8)-M(step1,5);
    dy = M(step1,9)-M(step1,6);
    dz = M(step1,10)-M(step1,7);
    dl = norm2([dx,dy,dz]);

    cheat.a(i) = atan2(dz,dy) *180/pi; %alpha (degrees)
    cheat.b(i) = asin(dx/dl) *180/pi;  %beta (degrees)
    
    D.x2{i} = x2;
    D.dE{i} = dEcur;

    clear x2
end

cheat.bcontained = bcontained;
cheat.pcontained = pcontained;
cheat.pEescaped = pEescaped;

D.img = img;
D.cheat = cheat;

if plotflag
    surf(img,'linestyle','none')
    view(2)
    axis equal
end



return  %end of program
%%


function x = SampleGauss2(s)
%x = SampleGauss2([N,M])
r = rand(s(1),s(2));
x = sqrt(2)*erfinv(r.*2-1);

function out = norm2(x)
%function out = norm2(x)
% x is three columns.
% return vector normal of each 3-vec.
out = sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2);
