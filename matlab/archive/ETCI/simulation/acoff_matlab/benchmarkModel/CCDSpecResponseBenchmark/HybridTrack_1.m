function Track = HybridTrack(img,lt,dEref,cheat,plotflag)
%function Track = HybridTrack(img,lt,dEref,cheat,plotflag)
%
%Algorithm for finding initial electron trajectory.
%Hybrid of thinning method (bwmorph) and stepping method.
%
%Inputs:
% img: matrix of energy values. (electron track image)
% lt: low threshold for making binary image
% dEref: table of dE/dx_ref values (col1: energy, col2: dE/dx_ref)
% plotflag: 1 = plot result.
% cheat: structure of info from DiffuseTrack. See below in Outputs.
%
%Output:
% Track: a structure of information on the track:
%  Track.img: copy of input image
%  Track.lt: low threshold input, keV.
%  Track.Etot: total energy of track = sum(img(:)); in keV.
%  Track.thin: binary image of thinned track
%  Track.ends: number of ends found on thinned image
%  Track.Eend: energy measured at end; used to determine starting end.
%  Track.x and Track.y: coordinates of points along an end of the track.
%       1 is the start of the track. In units of pixels.
%  Track.w: FWHM of track at each point, in um. (corresponds to x and y.)
%  Track.dE: dE/dx at each point (sum over cut), in keV/um. (corresponds to x and y.)
%  Track.a0: alpha at each point, units of degrees. (corresponds to x and y.)
%  Track.alpha: alpha angle measured for this track, units of degrees.
%       (azimuthal component of the initial electron trajectory)
%  Track.cosbeta: cosine of beta angle measured for this track.
%       cosbeta = dEref / median(dE(1:n)
%       (polar angle, + or - from pixel plane)
%  Track.ThSS: theta step size; defines units of a0 and alpha.
%
%  Track.cheat: additional information given from DiffuseTrack simulation
%       .cheat.Etot: actual electron energy (keV)
%       .cheat.Edep: actual energy deposited by electron (before noise) (keV)
%       .cheat.contained: logical for track contained in detector
%       .cheat.a: true value of alpha (degrees)
%       .cheat.b: true value of beta (degrees)
%       .cheat.z: starting depth from pixel plane (um)
%       .cheat.x: x-coordinate of starting pixel (pixels)
%       .cheat.y: y-coordinate of starting pixel (pixels) Track.img(x,y)



%The lowthresh should cut out some diffusion,
% so that the thinned track is closer to the actual track,
% and the correct end can be identified more accurately.
%If the lowthresh is too high, the initial end may be cut off.

imgbw = logical(img>lt);    %create binary image
% img(~imgbw) = 0;           %ignore small energy depositions under the low threshold (for model tracks)
imgbw = imgbw + 0;          %convert to numerical array
thin = bwmorph(imgbw,'thin',inf);   %thin image
conn = ones(3);                     %connectivity matrix 111/111/111
   
nb = conv2(single(thin),conn)-1;            %2D convolution to count neighbors, less one
nb = nb(2:size(nb,1)-1,2:size(nb,2)-1);     %crop border after convolution expands image
nb = nb .* thin;                           %only consider pixels in thin track

e = nb .* (nb==1);                  %pick out ends

%find ends
fe = find(e);
Track.ends = length(fe);
if length(fe) == 1
    Track.err = '1 end found';
    return
elseif isempty(fe)
    Track.err = 'no ends found';
    %try a different low threshold...
    %
    return
end

%Trying to identify starting end
%1.0
%Sum energy from pixels with R=2 pixels of end pixel.
convconn = [0,1,1,1,0; 1,1,1,1,1; 1,1,1,1,1; 1,1,1,1,1; 0,1,1,1,0];
Eend = conv2(img,convconn);    %sums energy over c
Eend = Eend(3:size(Eend,1)-2,3:size(Eend,2)-2);
Eend = Eend(logical(e)); 
%sort from min to max
[Eend,imin] = sort(Eend);   %new Eend = old Eend(imin)


%initialize
is = 1; intersection = zeros(length(imin),2);
xg=0; yg=0;
thin2 = thin;

%% Determine starting end
for i=1:length(imin)
    %check to see if this end is "good":
    %check length to intersection in thin. if it is too close to an intersection, the end
    %  may just be a diffuse elbow that was narrowed into an end by 'thin'. this effect
    %  will vary with the low threshold.
    [xo(i,1),yo(i,1)] = ind2sub(size(e),fe(imin(i)));   %coordinates of end
    x=xo(i,1); y=yo(i,1);
%     nconn = [1,1,1; 1,0,1; 1,1,1];
%     curend = zeros(size(img),'int8');
%     curend(x,y) = 1;
%     n = conv2(curend,nconn).*thin;
%     n = n(2:size(n,1)-1,2:size(n,2)-1);
    
    
        
    for j=2:5   %points along end of thin
        thin2(x,y) = 0;
        if sum([thin(x-1,y-1);thin(x,y-1);thin(x+1,y-1);thin(x+1,y);...
                thin(x+1,y+1);thin(x,y+1);thin(x-1,y+1);thin(x-1,y)]) > 2
            %this is an intersection
            %.....
            fi = find(intersection(1:is-1,1)==x & intersection(1:is-1,2)==y);
            if ~isempty(fi)
                %ignore shorter end
                if j < intersection(fi,3)
                    %use this end
                    xg = xo(i,:); yg = yo(i,:);
                else
                    %use other end
                    xg = xo(fi,:); yg = yo(fi,:);
                end
            else
                intersection(is,1:3) = [x,y,j];
            end
            break
        else
            if thin2(x-1,y-1)
                x=x-1; y=y-1;
            elseif thin2(x,y-1)
                y=y-1;
            elseif thin2(x+1,y-1)
                x=x+1; y=y-1;
            elseif thin2(x+1,y)
                x=x+1;
            elseif thin2(x+1,y+1)
                x=x+1; y=y+1;
            elseif thin2(x,y+1)
                y=y+1;
            elseif thin2(x-1,y+1)
                x=x-1; y=y+1;
            elseif thin2(x-1,y)
                x=x-1;
            end
        end
        xo(i,j)=x; yo(i,j)=y;
    end
    
    if xg~=0
        break
    end
    
    if intersection(is,1)==0
        %this track end is good
        xg = xo(i,:); yg = yo(i,:);
        break
    else
        is = is+1;
    end
    
end

if plotflag
    figure;
    subplot(2,2,1); surf(img,'linestyle','none'); view(0,90); axis equal; colormap('hot')
    subplot(2,2,2); surf(imgbw,'linestyle','none'); view(0,90); axis equal;
    subplot(2,2,3) 
    surf(thin+0,'linestyle','none'); view(0,90)
    hold on
    plot3(yg'+0.5,xg'+0.5,50.*ones(length(xg),1),'.g')
    axis equal
    pause(2)
    hold off
    subplot(2,2,4)
end

clear x y xo yo is 

%% Switch over into stepping mode.
%  now that we have [xg',yg'] coordinates of first 5 pixels of end of thinned track.

%Pixel size: to calculate dE/dx in keV/um
pixsize = 10.5; %microns
%Position step size
PSS = 2.625 / pixsize;     %units are pixels
%Cut step size
CSS = 2.625 / pixsize;     %pixels
%Cut (total) width
CW = 105 / pixsize;         %pixels

%Cut angle step size - should be factor of 45
ThSSd = 3;   %ThSSd in degrees
%Search angle - should be multiple of 2
%  max angle change in 1 step = SA/2
SA = 48/ThSSd;  %units are ThSS
%Work in units of ThSS instead of radians or degrees:
maxind = 360/ThSSd; pi0 = maxind/2; %pi in units of ThSS
ThSS = ThSSd*pi/180;     %ThSS in radians
th(1:maxind) = ThSS:ThSS:maxind*ThSS;
thd(1:maxind) = ThSSd:ThSSd:maxind*ThSSd;
% th(i) converts angle i (an integer) from ThSS to radians.
% thd(i) converts angle i (an integer) from ThSS to degrees.

%Low threshold for end of track
LT = 0.1;   %keV
%Cut low threshhold: for truncating cut width   [not implemented yet]
CLT = 0.5;  %keV
%Interpolation method for cuts
method = 'linear';

%Starting coordinates
x(1) = xg(length(xg));
y(1) = yg(length(yg));
%Cut for rotation
cut0y = -CW/2:CSS:CW/2;
cut0x = zeros(size(cut0y));     %Cut is along y-axis so that is perpendicular to trajectory
%distance from center to a point
d = abs(cut0y);
%Pre-rotated cuts
cut0 = cell(1,2*pi0);
for i=1:2*pi0
    R = [cos(th(i)),sin(th(i)); -sin(th(i)),cos(th(i))];
    cut0{i} = [cut0x',cut0y']*R;    %CCW rotation.
end



%% Find the first step angle. th is direction of step; cut is perpendicular
%look at two pixels of thin to guess the initial th0 direction
dx = xg(length(xg)-1)-xg(length(xg));
dy = yg(length(yg)-1)-yg(length(yg));
if dx==0 & dy==0
    Track.err = 'First Step angle failure';
    return;
end

th0 = atan(dy/dx) / ThSS;
%th0 will be a multiple of 45 degrees, thus integer in ThSS
if dx < 0
    th0 = th0 + pi0;
end

%Check widths
thc = th0-SA/2:th0+SA/2;    %length = SA+1
%place all angles into [0,2*pi0]:
thc(logical(thc>2*pi0)) = thc(logical(thc>2*pi0)) - 2*pi0;
thc(logical(thc<1)) = thc(logical(thc<1)) + 2*pi0;
cutE = cell(1,SA+1); wm = zeros(1,SA+1); %initialize
for j=1:SA+1
    xcur = x(1) + cut0{thc(j)}(:,1);
    ycur = y(1) + cut0{thc(j)}(:,2);
    %logcur is a logical which excludes out-of-bounds points.
    logcur = ~logical(xcur > size(img,1) | xcur < 1 | ycur > size(img,2) | ycur < 1);
    xcur = xcur(logcur);
    ycur = ycur(logcur);
    cutE{j} = interp2(img,ycur,xcur,method);   %interp2 switches x,y ...
    cutE{j}(logical(cutE{j}<0)) = 0;
    %--insert cut low thresh code here--
    wm(j) = sum(cutE{j} .* d(logcur)');     %(energy*distance) metric
end
[asdf,ind] = min(wm);   %get index of minimum value of width metric
th0(1) = thc(ind);  %units of ThSS

%Fit gaussian width and save FWHM
wfit = fit(cut0y',cutE{ind},'gauss1',...
    'StartPoint', [max(cutE{ind}),0,CW/5], ...
    'Lower', [0.5*max(cutE{ind}),-CW/2,.25], 'Upper', [2*max(cutE{ind}),CW/2,CW/2]);
w(1) = 1.665 * wfit.c1 * pixsize;       %FWHM in um
dE(1) = sum(cutE{ind}) * CSS / pixsize; %not really dE, but dE/dx in keV/um.

%Adjust x,y to centroid
x(1) = sum(cutE{ind}.*(x(1)+cut0{thc(ind)}(:,1))) / sum(cutE{ind});
y(1) = sum(cutE{ind}.*(y(1)+cut0{thc(ind)}(:,2))) / sum(cutE{ind});

if plotflag
    surf(img,'linestyle','none'); view(0,90)
    hold on
    plot3(y(1)+.5,x(1)+.5,50,'.c')
    axis equal
end

i=1;
%% Iterate.
while(1)
    i=i+1;
    %Take a step
    x(i) = x(i-1) + PSS*cos(th0(i-1)*ThSS);
    y(i) = y(i-1) + PSS*sin(th0(i-1)*ThSS);
    %Check low threshold
    if interp2(img,y(i),x(i),'nearest') < LT   %interp2 switches x,y ...
        x = x(1:i-1);
        y = y(1:i-1);
        break
    end
    
    %Check widths
    thc = th0(i-1)-SA/2:th0(i-1)+SA/2;
    thc(logical(thc>2*pi0)) = thc(logical(thc>2*pi0)) - 2*pi0;
    thc(logical(thc<1)) = thc(logical(thc<1)) + 2*pi0;
%     cutE = cell(1,length(th)); wm = zeros(1,length(th));
    for j=1:SA+1
        xcur = x(i) + cut0{thc(j)}(:,1);
        ycur = y(i) + cut0{thc(j)}(:,2);
        %logcur is a logical which excludes out-of-bounds points.
        logcur = ~logical(xcur > size(img,1) | xcur < 1 | ycur > size(img,2) | ycur < 1);
        xcur = xcur(logcur);
        ycur = ycur(logcur);
        cutE{j} = interp2(img,ycur,xcur,method);   %interp2 switches x,y ...
        cutE{j}(logical(cutE{j}<0)) = 0;
        %--insert cut low thresh code here--
        wm(j) = sum(cutE{j} .* d(logcur)');     %(energy*distance) metric
    end
    [asdf,ind] = min(wm);   %get index of minimum value of width metric
    th0(i) = thc(ind);

    %Fit gaussian width and save FWHM
    wfit = fit(cut0y',cutE{ind},'gauss1',...
        'StartPoint', [max(cutE{ind}),0,CW/5], ...
        'Lower', [0.5*max(cutE{ind}),-CW/2,.25], 'Upper', [2*max(cutE{ind}),CW/2,CW/2]);
    w(i) = 1.665 * wfit.c1 * pixsize;       %FWHM in um
    dE(i) = sum(cutE{ind}) * CSS / pixsize; %not really dE, but dE/dx in keV/um.
    
    %Adjust x,y to centroid
    x(i) = sum(cutE{ind}.*(x(i)+cut0{thc(ind)}(:,1))) / sum(cutE{ind});
    y(i) = sum(cutE{ind}.*(y(i)+cut0{thc(ind)}(:,2))) / sum(cutE{ind});

    %Check for an infinite loop in the track
    %If we are within (PSS/2) of a previous point, there is a problem.
    if min((x(i)-x(1:i-1)).^2 + (y(i)-y(1:i-1)).^2) < PSS^2/4
        Track.err = 'Infinite loop';
        return
    end
        
    if plotflag
        plot3(y(i)+.5,x(i)+.5,50,'.c')
        axis equal
    end

end

%% Calculations on the track

%reverse indices of all point-by-point variables
%  so that 1 is the start of the track.
x = x(i-1:-1:1);
y = y(i-1:-1:1);
w = w(i-1:-1:1);
th0 = th0(i-1:-1:1) - pi0;
dE = dE(end:-1:1);

Etot = sum(img(:));

%Measure width of track so we know how many points to skip.
wimax = min([8,length(w)]);
w1 = mean(w(1:wimax));
s = (w1 - pixsize) * 4/pixsize; %logbook 11/1/09

%calculate x prime and y prime for this energy
xp = sqrt(0.0825*Etot + 15.814) - 3.4;
xp = max([xp,0]); %do not let x prime < 0
yp = 2*xp + 3.4;
yp = max([yp,0]);

%Zeroth estimate of beta: guess 45.
beta0 = 45 *pi/180;

x0 = ceil(xp * cos(beta0) + s);
x0 = min([x0,length(dE)]);  %do not let parameters overshoot vector length
y0 = ceil(yp * cos(beta0) + s);
y0 = min([y0,length(dE)]);

%First estimate of beta.
cosbeta1 = interp1(dEref(:,1), dEref(:,2), Etot) / median(dE(x0:y0));
if cosbeta1 < 1
    beta1 = acos(cosbeta1);
else
    beta1 = 0;
end

x1 = ceil(xp * cos(beta1) + s);
x1 = min([x1,length(dE)]);
y1 = ceil(yp * cos(beta1) + s);
y1 = min([y1,length(dE)]);

%Second estimate of beta.
cosbeta2 = interp1(dEref(:,1), dEref(:,2), Etot) / median(dE(x1:y1));
if cosbeta2 < 1
    beta2 = acos(cosbeta2);
else
    beta2 = 0;
end

%Estimate of alpha.
alpha = median(th0(x1:y1))*ThSSd;   %in degrees

%Cheat...
if ~isempty(cheat)
    cosbeta1c = interp1(dEref(:,1), dEref(:,2), cheat.Etot) / median(dE(x0:y0));
    if cosbeta1c < 1
        beta1c = acos(cosbeta1c);
    else
        beta1c = 0;
    end
    x1c = ceil(xp * cos(beta1c) + s);
    x1c = min([x1c,length(dE)]);
    y1c = ceil(yp * cos(beta1c) + s);
    y1c = min([y1c,length(dE)]);

    cosbeta2c = interp1(dEref(:,1), dEref(:,2), cheat.Etot) / median(dE(x1c:y1c));
    if cosbeta2c < 1
        beta2c = acos(cosbeta2c);
    else
        beta2c = 0;
    end

    alphac = median(th0(x1c:y1c))*ThSSd; %in degrees
end

if plotflag; hold off; end

%% Format & save.

Track.img = img;
Track.thin = thin;
Track.x = x;
Track.y = y;
Track.w = w;
Track.a0 = (th0).*ThSSd;      %degrees
Track.alpha = alpha;    %already in degrees
Track.beta = beta2 *180/pi; %degrees
if ~isempty(cheat)
    Track.alphac = alphac;
    Track.betac = beta2c *180/pi;
end
Track.Etot = Etot;
Track.ThSS = ThSS;
Track.dE = dE;

%for optimizations...
Track.ends = Track.ends;
Track.lt = lt;
Track.Eend = Eend;

%% Cheat!

if ~isempty(cheat)
    Track.cheat = cheat;
    %already added to cheat structure:
    %cosbeta using actual electron energy (e.g. if detector captured full
    %energy of electron)

    Track.cheat.enddist = sqrt((Track.x(1)-cheat.x)^2 + (Track.y(1)-cheat.y)^2);  %2D distance from first algorithm point to true electron start

end