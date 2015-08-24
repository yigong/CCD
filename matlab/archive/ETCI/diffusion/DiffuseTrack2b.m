function [D] = DiffuseTrack2(M,psft,th,plotflag)
%function [D] = DiffuseTrack2(M,psft,th,plotflag)
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
% Coordinate system: defined by driftdim below (pixel face = 0, back face = (th/1000))
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
%D.img:     data matrix of ccd plane pixels.    (keV)
%   ***if segment_flag is used in settings below, D.img is a CELL ARRAY instead of a matrix.
%       each cell is one track matrix segmented by the segment code.
%D.cheat:
%  cheat.Etot:    actual electron energy        (keV)
%  cheat.Edep:    energy deposited in detector  (keV)
%  cheat.a:       alpha                         (degrees)
%  cheat.b:       beta                          (degrees)
%  cheat.z:       zp (starting z)               (um)
%  cheat.contained: logical for all energy deposited in active volume
%               (due to rounding errors, Edep ~= Etot)
%  cheat.x:     x coordinate of start of trajectory (pixels)
%  cheat.y:     y coordinate of start of trajectory (pixels)
%D.x2:      list of 3D coordinates that have been diffused. z=depth.
%D.dE:      list of energies diffused. corresponds to D.x2.

pixsize = 10.5;     %um
buffer = 8.5*pixsize;    %um            %point (0,0) is in the center of a pixel

N = max(M(:,1));    %number of particles

%which dimension is the drift of charge carriers (perpendicular to pixels)?
driftdim = 'x';

%apply segmentation code? (CCDsegment3.m)
segment_flag = false;

%apply noise? (10/13/2011)
noiseflag = false;

%% Particle info.

f1 = zeros(1,N);
% f2 = zeros(1,N);
q = zeros(1,N);
Ep = zeros(1,N);
p = zeros(1,N);

asdf=0; %dummy variable for Matlab 2008b

for i=1:N
    f1(i) = find(M(:,1)==i,1,'first');  %first record of this particle
%     f2(i) = find(M(:,1)==i,1,'last');   %last record of this particle
    q(i) = M(f1(i),4);          %charge
    Ep(i) = M(f1(i),14);        %initial energy
    p(i) = M(f1(i),2);          %parent ID
end

%% Remove multiple Compton interactions.
% Only look at first Compton electron.
fCompton = find(p==1);  %expect #1 to be the photon...
if length(fCompton)>1
    [asdf,mC] = min(M(f1,1));
    fCompton(mC) = [];  %remove the first Compton electron from this list.
    
    %remove the others from the matrix.
    i=1;    %index for fCompton list
    while i <= length(fCompton)
        %flag daughters of this particle for removal
        fadd = find(p==fCompton(i));
        fCompton = [fCompton, fadd];
        %remove this particle
        M(M(:,1)==fCompton(i),:) = [];
        i = i+1;
    end
    
    %redefine some variables
    N = max(M(:,1));
    f1 = NaN(1,N);
    for i=1:N
        if ~isempty(find(M(:,1)==i,1,'first'))
            f1(i) = find(M(:,1)==i,1,'first');
        end
    end
end

%% Step info.

%Inter-record steps.
x = [0,0,0; (M(1:end-1,8:10) + M(2:end,8:10))./2];   %average of final positions from one record to next.
ss = [0;norm2(diff(M(:,8:10)))];                %step sizes from one record to the next.
dE = [0;-diff(M(:,13))];        %difference in final energies of each record.
%A given entry contains the difference from the previous line, to this line.

%Intra-record steps (at beginning of each particle track)
x(f1(~isnan(f1)),:) = (M(f1(~isnan(f1)),5:7) + M(f1(~isnan(f1)),8:10))./2;
ss(f1(~isnan(f1))) = norm2(M(f1(~isnan(f1)),5:7) - M(f1(~isnan(f1)),8:10));
dE(f1(~isnan(f1))) = M(f1(~isnan(f1)),14) - M(f1(~isnan(f1)),13);
%note: I have used [14]-[13] instead of [15] because I am expecting to see the energy of the secondary in its parent.

dE(M(:,4)==0) = M(M(:,4)==0,15);    %photons: dE represents binding energy of Compton-/photo-electron, 
x(M(:,4)==0,:) = M(M(:,4)==0,8:10);   %  which is deposited by x-ray very near the atom.

%% Deal with secondaries.
%must subtract energy carried by secondary particle, from energy lost by parent.

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

%% Unit conversion, coordinate transform.

x = x .*1e3;        %mm to um
dE = dE .*1e-3;     %eV to keV

if strcmp(driftdim,'x')
    %  x(:,1) is the depth dimension. ('x' => 'z')
    %  x(:,2) is the first pixel dimension. ('y' => 'x')
    %  x(:,3) is the second pixel dimension. ('z' => 'y')
    x2(:,3) = x(:,1);
    x2(:,1) = x(:,2);
    x2(:,2) = x(:,3);
    x = x2; clear x2
elseif strcmp(driftdim,'z')
    % no action needed.
else
    error('Drift dimension not defined well... y not supported yet')
end

escapes = x(:,3) > th | x(:,3) < 0;

if any(escapes & M(:,4)==-1)
    contained = false;
else
    contained = true;
end

x = x(~escapes & dE>0,:);     %without dE>0, empty energy deposition (e.g. from movement of source photon) might cause the track image to be huge
dE = dE(~escapes & dE>0);

%% Remove stray points
%  (e.g. distant brehmsstrahlung deposition, ???)

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
    x_0 = x2(1,1) / 10.5;
    y_0 = x2(1,2) / 10.5;
    
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

%% Set up image
        
% zres = 0.5; %um resolution of PSFT
zres = th/(length(psft)-1);

int=0.5; %um
%This is the resolution of the fine grid. Energy deposition coordinates
%will be rounded off to a grid position.

%int must match the xres of PSFT

%initialize fine grid of track image.
fgrid=zeros(imsize.*pixsize./int);

gs = size(psft{1},1);    %gaussian size of psf.
gs2 = (gs-1)/2;

%lookup psf for each energy deposition point.
%sum into image of resolution 'int'.
for i=1:size(x2,1)
    zi = round(x2(i,3)/zres)+1;     %index of psft to use
    Espread = dE(i,1).*psft{zi};    %(already keV units)
    xi = round(x2(i,1)/int);        %round x,y to the nearest fine grid coordinate
    yi = round(x2(i,2)/int);
    fgrid(xi-gs2:xi+gs2,yi-gs2:yi+gs2) = fgrid(xi-gs2:xi+gs2,yi-gs2:yi+gs2) + Espread;
end

%create CCD image at pixel resolution.
img=zeros(imsize);
for p=1:imsize(1)
    for r=1:imsize(2)
%         SGrid(p,r)=sum(sum(Grid(51+pxsize*(p-1)/int:50+pxsize*(p)/int,51+pysize*(r-1)/int:50+pysize*(r)/int)));  
        img(p,r) = sum(sum(fgrid(1+(p-1)*pixsize/int:p*pixsize/int,1+(r-1)*pixsize/int:r*pixsize/int)));
    end
end

if noiseflag
    %Noise: dark current + Fano
    F = 0.14;           %Fano factor to use. Silicon: 0.08 < F < 0.14 (Knoll)
    blsig = 0.01874;        %std dev of background level, in keV. Nov 2007 data: 0.01874
    epcc = 3.7e-3;      %energy per eh pair, in keV
    cc = img ./ epcc; %average # charge carriers in each pixel
    u = sqrt(blsig.^2 + epcc^2 .* F .* cc);    %add stdev's in quadrature to get total sigma of noise

    img = img + u.*SampleGauss2(size(img));     %add gaussian noise of magnitude u over image
end

step1 = f1(find(q==-1,1));
%Set up outputs
cheat.Edep = sum(dE);     %already in keV
cheat.Etot = M(step1,14) .* 1e-3; %first E_initial, in keV
%due to rounding errors, abs(Etot - Edep) = ~4e-7 *Etot, for contained tracks
cheat.contained = contained;
cheat.x = x_0;  %in pixel coordinates of image
cheat.y = y_0;  %in pixel coordinates of image
cheat.z = x2(1,3);  %already in um


dx = M(step1,8)-M(step1,5);
dy = M(step1,9)-M(step1,6);
dz = M(step1,10)-M(step1,7);
dl = sqrt(dx^2 + dy^2 + dz^2);

cheat.a = atan2(dz,dy) *180/pi; %alpha (degrees)
cheat.b = asin(dx/dl) *180/pi;  %beta (degrees)

D.img = img;
D.cheat = cheat;

D.x2 = x2;
D.dE = dE;

if plotflag
    surf(img,'linestyle','none')
    view(2)
    axis equal
end

%% Segmentation code.
if segment_flag
    %really would like to simulate x,y position too, for imaging...
    %but that's a bit beyond this version of the code, i think.
    
%     subpar = [];
%     [T,E] = CCDsegment3(img,0.55,'Subarray',subpar);
    [T,E] = CCDsegment3(img,0.55);  %default options
    for i=1:length(T)
        D.img{i} = T{i}.img;
    end
else
    D.img = img;
end


function x = SampleGauss2(s)
%x = SampleGauss2([N,M])
r = rand(s(1),s(2));
x = sqrt(2)*erfinv(r.*2-1);

function out = norm2(x)
%function out = norm2(x)
% x is three columns.
% return vector normal of each 3-vec.
out = sqrt(x(:,1).^2 + x(:,2).^2 + x(:,3).^2);
