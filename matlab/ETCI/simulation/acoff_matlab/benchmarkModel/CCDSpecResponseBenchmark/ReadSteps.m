function out = ReadSteps(M,x1,x2)
%function out = ReadSteps(M,x1,x2)
%
%M = geant4 output data for photon scattering
%x1 = minimum x-value of silicon in geometry (e.g. -21.5)
%x2 = maximum x=value of silicon in geometry (e.g. 10.5)

%{
columns:
 1      track ID
 2      step #
 3      parent ID
 4      charge
 5,6,7  init pos
 8,9,10 final pos
11      accumulated track length
12      step length
13      E final energy
14      dE change in energy
15      E+dE initial energy
16      Edep energy deposited

%}

%.mat files are faster to load
%{
try
    load([filename,'.mat']);
catch
    M = load(filename);
    save([filename,'.mat'],'M');
end
%}
n = size(M,1);

%% USEFUL:
len = size(M,1);

%step #1
fs1 = logical(M(:,2)==1);

%interactions = energy deposited
fed = logical(M(:,14)~=0);  %uses DE

%full energy depositions
ffed = logical(M(:,13)==0); %E_final

%interactions in silicon
% % x1 = -21.5; x2 = 10.5;
y1 = -17.5; y2 = 17.5;
z1 = y1; z2 = y2;
fsi = logical(M(:,8)>=x1 & M(:,8)<=x2 & M(:,9)>y1 & M(:,9)<y2 & M(:,10)>z1 & M(:,10)<z2);
fsi = fsi & fed;

%interactions in CZT
fcz = logical(M(:,8)>x2 | M(:,9)<y1 | M(:,9)>y2 | M(:,10)<z1 | M(:,10)>z2);
fcz = fcz & fed;

%first interactions
ffi = fed(2:len) & fs1(1:len-1);
ffi = logical([0;ffi]);

%last interactions
fli = fed(1:len) & logical([fs1(2:len);1]);

%first interaction in silicon
ffisi = ffi & fsi;

%polar angle of first interaction
th = atan(sqrt((M(ffi,9)-M(ffi,6)).^2 + (M(ffi,10)-M(ffi,7)).^2)./(M(ffi,8)-M(ffi,5)));
th(logical(th<0)) = th(logical(th<0)) + pi;

%count # interactions
nint = fed+0;
nint(:,2) = fed .* (1 + [0;nint(1:len-1,1)]);
while(~all(nint(:,2)==nint(:,1)))
    nint(:,1) = nint(:,2);
    nint(:,2) = fed .* (1 + [0;nint(1:len-1,1)]);
end
nint = nint(:,2);
%nint = nint(fli);      %the entries in nint are the number of interactions for each photon


% ffisi & nint
t1 = ffisi;     %temp
n2 = zeros(len,1);
t1(logical(ffisi & fli)) = 0;    %kill if only one interaction
while(any(t1))
    t1 = logical([0;t1(1:len-1)]);   %shift down one
    n2 = nint .* logical(t1 & fli | n2>0);   %save
    t1(logical(t1 & fli)) = 0;      %kill at end
end
fmssi = n2 .* nint;

%ffisi & nint & ffed
fmssife = fmssi .* ffed;

%count # interactions in Si
%count # interactions in CZT

%total detector energy spectrum

%% OUTPUT
N = sum(fs1+0);
fullE = sum(logical(fmssife)+0);
int2 = sum(logical(fmssi)+0);
intany = sum(ffi+0);


out.N = N;
out.intany = intany;
out.int2 = int2;
out.fullE = fullE;
% 
% disp(['Total number of photons: ',num2str(N)])
% disp(['Photons with some interaction: ',num2str(intany),' +- ',num2str(sqrt(intany))])
% disp(['Photons with >1 interaction; first interaction in Si: ',num2str(int2),' +- ',num2str(sqrt(int2))])
% disp(['Photons with >1 interaction; first interaction in Si; full E deposited: ',num2str(fullE),' +- ',num2str(sqrt(fullE))])

%% TROUBLESHOOTING:
%{
dx = M(:,8) - M(:,5);
dy = M(:,9) - M(:,6);
dz = M(:,10) - M(:,7);

%first interaction associated with backwards-directed vector
ft1 = ffisi & logical(dx<0);
%}
return

%{
%% OLD STUFF
%{
%.mat files are faster to load
try
    load([filename,'.mat']);
catch
    M = load(filename);
    save([filename,'.mat'],'M');
end

n = size(M,1);
%}

%find first interactions in silicon
%from pos0:
fsi2 = find(M(2:n,5)>=-21.5 & M(2:n,5)<11.5 & ... %front to back (x) coordinates of silicon
    M(2:n,6)<17.5 & M(2:n,6)>-17.5 & ...        %y coordinates
    M(2:n,7)<17.5 & M(2:n,7)>-17.5 & ...        %z coordinates
    M(1:n-1,5)==-1016.5);                        %previous line is 0 energy deposited
                                                %  means that this is the first interaction

                                                
%find first interactions in silicon
%from pos1:
fsi = find(M(2:n,8)>-21.5 & M(2:n,8)<=11.5 & ... %front to back (x) coordinates of silicon
    M(2:n,9)<=17.5 & M(2:n,9)>=-17.5 & ...        %y coordinates
    M(2:n,10)<=17.5 & M(2:n,10)>=-17.5 & ...        %z coordinates
    M(1:n-1,5)==-1016.5);                        %previous line is 0 energy deposited

%find first interactions in silicon
%from pos_av:
fsi = find(M(2:n,8)>-21.5 & M(2:n,8)<=11.5 & ... %front to back (x) coordinates of silicon
    M(2:n,9)<=17.5 & M(2:n,9)>=-17.5 & ...        %y coordinates
    M(2:n,10)<=17.5 & M(2:n,10)>=-17.5 & ...        %z coordinates
    M(1:n-1,5)==-1016.5);                        %previous line is 0 energy deposited


%full energy depositions
ffed = find(M(1:n,13)==0);      %E_final = 0
                                                    
%no interaction
fno = find(M(2:n,5)==-1016.5 & M(1:n-1,5)==-1016.5);
% a line of no energy deposited only means no interaction if there is not an energy
% deposit on the next line

%some interaction
fso = find(M(2:n,5)==-1016.5 & M(1:n-1,14)~=0);
% energy deposited on one line but a new photon on the next

%}