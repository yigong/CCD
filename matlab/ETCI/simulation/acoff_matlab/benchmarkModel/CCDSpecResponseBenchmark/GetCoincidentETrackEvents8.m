function Res = GetCoincidentETrackEvents8(filename)
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

    % Start Condition -- (This will need updating for Radioacitve Decay
    % Sources 
    % Photon Source
        startX = S(:,1)+S(:,2)+S(:,3);
    % Index Location of Start Steps
        FE = find(startX==2);
    
    % Radioactive Sources: (Start of large -x position)....
        % startX = ...
        % FE = find(startX < -10);  % starting Source Postion, x = - 10 cm
    
    NRead = length(FE); %Find total number of Events
    disp(['Number of Read Events: ' num2str(NRead)]);
    
% Set Containters 
ncoin = 0;
nfound = 0;
S_old = [];

for i = 1:NRead
    
   NE = [FE(2:end);size(S,1)]-FE;  %Find Number of entrys per event
   FE2 = FE(find(NE > 2));         %Filter events that are only single interactions
   NE2 = NE(find(NE > 2));         %Using NE2 for offsetting...
    %Loop over FE2 to find  Events of One's Choosing.
    % 10/18/12 - No Event Selection Currently
    % If want to add Event Selections qualify F1, F2... accordingly, and
    % re-add 'if ~isempty(Fn)' statement
    
        S2 = S(FE(i):FE2(i)+NE2(i)-1,:);
        %F1 = find(S2(:,4) == 0 & S2(:,7) > -10 & S2(:,15)>0); % Track in si-CCD
        %F2 = find(S2(:,4) == 0 & S2(:,7) < -10 & S2(:,15)>0); % Track in Ge
        %if ~isempty(F1) %& ~isempty(F2) % Choose if want si-CCD only tracks or Coinc. Tracks 
            nfound = nfound + 1;
            ncoin = ncoin + 1;
            Res{ncoin} = S2;
        %end   
    if (mod(i,100) == 0)
        disp(['Found ',num2str(nfound),' events, Total events = ',num2str(ncoin)]);
        nfound = 0;
    end
end
disp(['Found ',num2str(nfound),' events, Total events = ',num2str(ncoin)]);

% Save Track Matrix, (optional)
% F = strfind(filename,'.dat');
% sname = [filename(1:F(1)-1),'_c.mat'];
% save(sname,'Res');
    
    
    