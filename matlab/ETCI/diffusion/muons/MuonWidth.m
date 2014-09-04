function out = MuonWidth(Track,MuonCoord,MuonWidth,plotflag)
%function out = MuonWidth(Track,MuonCoord,MuonWidth,plotflag)
%
% Analyze widths of muon tracks (experimental data). Ignore spurs, crossed tracks, etc.
%
%Inputs:
%
% Track         cell array of track structures
% MuonCoord     m by 4 array of coordinates. Each row i gives [x1 y1 x2 y2] for Track{i}
% plotflag      1 = plot widths; 0 = don't
% MuonWidth     cell array of previous result - so you don't have to re-compute
%
%Output:
%
% ...

if iscell(MuonWidth)
    st = 1 + length(MuonWidth);
    out = MuonWidth;
else
    st = 1;
end

z = min(length(Track),size(MuonCoord,1));

for i = st:z
    disp(num2str(i))
    if z>st
%         progressbar((i-st)/(z-st));
        
    end
    if ~isempty(Track{i}) && MuonCoord(i,1)==-1
        %flip track around
        disp(['i = ',num2str(i),' of ',num2str(z)])
        out{i} = OneMuonWidth(Track{i}.TrackImg,MuonCoord(i,4),MuonCoord(i,5),MuonCoord(i,2),MuonCoord(i,3),0);
        
        if plotflag
            %plot into widths window
            plot(out{i}.z,out{i}.w,'Color',[0,i/z,i/z]); hold on
        end
    elseif ~isempty(Track{i}) && MuonCoord(i,1)>0   %skip empty Track cells or 0s in coordinates
        disp(['i = ',num2str(i),' of ',num2str(z)])
        out{i} = OneMuonWidth(Track{i}.TrackImg,MuonCoord(i,2),MuonCoord(i,3),MuonCoord(i,4),MuonCoord(i,5),1);
        
        if plotflag
            %plot into widths window
            plot(out{i}.z,out{i}.w,'Color',[0,i/z,i/z]); hold on
        end
    end
end
% progressbar(1);

        


function out = OneMuonWidth(img,x1,y1,x2,y2,plotflag)
%function out = MuonWidth(img,x1,y1,x2,y2,plotflag)
%
% Analyze width of muon track (experimental data). Ignore spurs, crossed tracks, etc.
%
% img:      Track image. e.g. Track{k}.TrackImg
% x1,y1:    coordinates of narrow end of track.
% x2,y2:    coordinates of wide end of track.
%
% out.w:    widths in order from narrow to wide
% out.dE:   dE in order from narrow to wide
% out.z:    depths (linear from front to back of CCD) to correspond to w

maxdE = 5;
maxw = 6;
maxdw = .25;


%buffer image edges
buf = 5;
img2 = zeros(size(img,1)+2*buf,size(img,2)+2*buf);
img2(buf+1:size(img,1)+buf,buf+1:size(img,2)+buf) = img;
clear img
img = img2;
clear img2

%buffer input values
x1 = x1+buf;
y1 = y1+buf;
x2 = x2+buf;
y2 = y2+buf;


%Fit line y=mx+b
m = (y2-y1)/(x2-x1);
% b = y1 - m*x1;

%xy length
% len = sqrt((y2-y1)^2+(x2-x1)^2);

% %number of points
% N = ceil(len / 0.5);

if x1 < x2
    step = 0.5 / sqrt(1 + m^2);
else
    step = -0.5 / sqrt(1 + m^2);
end

%points along muon track
x = [x1:step:x2,x2];
y = [y1:m*step:y2,y2];

%slope of orthogonal cuts
mcut = -1 / m;

css = 0.25; %cut step size
cl = 7;     %cut length (total)
stepcut = css / sqrt(1+mcut^2);
xcut = -(cl/css/2)*stepcut : stepcut : (cl/css/2 + 0.5)*stepcut;    %0.5 to make sure rounding error doesn't cut it short
ycut = -(cl/css/2)*mcut*stepcut : mcut*stepcut : (cl/css/2 + 0.5)*mcut*stepcut;

% [xg,yg] = meshgrid(1:size(img,1),1:size(img,2));

% progressbar(0);
for i=1:length(x)
    cut(i,1:(cl/css+1)) = interp2(img,xcut+x(i),ycut+y(i),'cubic');
    cut(i,logical((cut(i,:)<0)|isnan(cut(i,:)))) = 0;
    xfit = -cl/2:css:cl/2;
    fitopts = fitoptions('Method','NonlinearLeastSquares','StartPoint',[max(cut(i)),0,2],'Lower',[0,1,0],'Upper',[1.5*max(cut(i))+1,cl/css+2,cl/css]);
    [f,gof] = fit(xfit',(cut(i,:))','Gauss1',fitopts);
    
    out.dE(i) = sum(cut(i,:)) * css;
    out.rmse(i) = gof.rmse;

    %reject based on width and dw
    if f.c1 <= maxw && (i==1 || ~(abs(f.c1 - out.w(i-1))>maxdw))   %ignore dw measure on first point - no w(i-1) exists
        out.w(i) = f.c1;
    else
        out.w(i) = NaN;
    end
    clear f
%     progressbar(i / length(x));
end

%reject based on dE - also reject neighbors
bad = logical(out.dE > maxdE)+0;
c = [1,1,1];
rej = conv(c,bad);
rej = logical(rej(2:length(rej)-1));
out.w(rej) = NaN;


out.z = (x - x1)./(x2-x1) .*650;    %in micron
out.img = img;

out.x1 = x1;
out.y1 = y1;
out.x2 = x2;
out.y2 = y2;

if plotflag
    surf(img,'linestyle','none');view(0,90)
    hold on
    plot3([x1;x2]+0.5,[y1;y2]+0.5,[100;100],'-r')
    for i=1:size(cut,1)
        if ~isnan(out.w(i))
            plot3(xcut'+x(i)+0.5,ycut'+y(i)+0.5,xcut'+150,'-r')
        end
    end
    axis equal

    hold off

end