function Res = GetCoincidentETrackEvents5(filename)


fname = filename;

%Dat File
Res = {};

if nargin == 0
    [fname, pname] = uigetfile( ...
       {'*.*','All Files (*.*)'}, ...
        'Pick a file');
else
    pname = [];
end

S = load(fname,'-ascii');

if S == -1
    error('Failed to open file');
    return;
end 

    %RAD = sqrt(S(:,5).^2 + S(:,6).^2 + S(:,7).^2);
    startX = S(:,1)+S(:,2)+S(:,3);
    %Start Condition
    FE = find(startX==2);
    %FE = find(RAD > 900);  %Old: Find start
    
    NRead = length(FE); %Find total number of Events
    disp(NRead);
    
% Set Containters 
ncoin = 0;
S_old = [];

for i = 1:NRead
    
   NE = [FE(2:end);size(S,1)]-FE;  %Find Number of entrys per event
   FE2 = FE(find(NE > 2));         %Filter events that are only single interactions
   NE2 = NE(find(NE > 2));         %Using NE2 for offsetting...
    %Loop over FE2 to find coincident events
    nfound = 0;
    
        S2 = S(FE(i):FE2(i)+NE2(i)-1,:);
        F1 = find(S2(:,4) == 0 & S2(:,7) > -10 & S2(:,15)>0); % Track in si-CCD
        F2 = find(S2(:,4) == 0 & S2(:,7) < -10 & S2(:,15)>0); % Track in Ge
        if ~isempty(F1) %& ~isempty(F2) % Choose if want si-CCD only tracks or Coinc. Tracks 
            nfound = nfound + 1;
            ncoin = ncoin + 1;
            Res{ncoin} = S2;
        end   
    %if (i . 10) == 0
        disp(['Found ',num2str(nfound),' events, Total events = ',num2str(ncoin)]);
    %end
end

% Save Track Matrix, (optional)
% F = strfind(filename,'.dat');
% sname = [filename(1:F(1)-1),'_c.mat'];
% save(sname,'Res');
    
    
    