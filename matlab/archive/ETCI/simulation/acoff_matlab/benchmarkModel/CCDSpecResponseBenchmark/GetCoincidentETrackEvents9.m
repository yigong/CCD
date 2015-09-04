function Res = GetCoincidentETrackEvents9(filename,radDecayFlag)
fname = filename;
if ~exist('radDecayFlag','var')
    radDecayFlag = false;
end
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
    
    if radDecayFlag
    % Radioactive Sources: (Start of large -x position)....
        startXYZ =  [S(1,5) S(1,6) S(1,7)];
        RAD0 = sqrt(S(1,5).^2 + S(1,6).^2 + S(1,7).^2);
        RAD = sqrt(S(:,5).^2 + S(:,6).^2 + S(:,7).^2);
        startX = (S(:,3)+S(:,4));
        activeCCDLayer = 2; % Column 17 is layer number
        CCDLayer1 = 1;
        CCDLayer2 = 2;
        CCDLayer3 = 3;
        electronCharge = -1;
        FE = find(RAD == RAD0);  % Find start of Photons
        %FE = find(startX == 0);  % starting Source Postion, x = - 10 cm
        %FE1 = find((S(:,4) == electronCharge)); % Only Look at Electrons
        FE1 = find((S(:,17) == activeCCDLayer)|(S(:,17) == CCDLayer1)|(S(:,17) == CCDLayer3));
        S_new = S(FE1,:);length(S_new)
        %clear S;
        % Start of an electron  (Don't want Electron 'finding based on
        % 1+-1== 0... assumes start of electron starts in active layer which is 99% true
        % Base seperation on change of charge from 0 to -1 
        %startX = (S_new(:,3)+S_new(:,4));
        %startElec = (S_new(:,3)+S_new(:,4));
        % Start of Photon - Photon hitting CCD: 
        charge = S_new(:,4); % Photon ==0, electron == -1
        length(charge)
        chg1 = charge(1:(length(charge)-1)); length(chg1)
        chg2 = charge(2:length(charge)); length(chg2)
        chargeChange = chg2-chg1; %chg1-chg2;
        length(chargeChange)
        chargeChange = [1;chargeChange];length(chargeChange)
        % Start of CCD Event
        %FE = find(startX == 0);
        %FE = find(charge == 0); %Assumes only 1 photon causes 1 
                                     %electron tracks 'bunch' (temp solution) 
                                     %(will miss photons from multiple compton scattering cases...)
        FE = find(chargeChange == 1);
        length(FE)
        %FE = [FE,0];
        %length(FE)
        
    else
        % Start Condition:
        % Sources 
        % Photon Source
            startX = S(:,1)+S(:,2)+S(:,3);
        % Index Location of Start Steps
            FE = find(startX==2);
    end
    
    NRead = length(FE); %Find total number of Events
    disp(['Number of Read Events: ' num2str(NRead)]);
    
% Set Containters 
ncoin = 0;
nfound = 0;
S_old = [];


    if radDecayFlag
        % Working on Sorting for Rad-Decay Sources
        %NE = [FE(2:end);size(S,1)]-FE;  %Find Number of entrys per event 
        %NE = [FE(2:end);size(S_new,1)]-FE; % NEW
        NE = [FE(2:end);size(S_new,1)]-FE; % NEW Newer
        FE2 = FE(find(NE >= 2));         %Filter events that are only single interactions (skips loner photons)
        NE2 = NE(find(NE >= 2));         %Using NE2 for offsetting...
        radNRead = length(NE2); %Find total number of Events
        disp(['Number of Electrons Read Events: ' num2str(radNRead)]);
        for j = 1:radNRead
            
            %S2 = S(FE(j):FE2(j)+NE2(j)-1,:); % old
            S2 = S_new(FE2(j):(FE2(j)+NE2(j)-1),:); % NEW 
            % Within FE2, find First Start of each electron...
            nfound = nfound + 1;
            ncoin = ncoin + 1;
            Res{ncoin} = S2;
                    if (mod(j,1000) == 0)
                        disp(['Found ',num2str(nfound),' events, Total events = ',num2str(ncoin)]);
                        nfound = 0;
                    end
        end 
    else
            % Working Sorting for Photons
            NE = [FE(2:end);size(S,1)]-FE;  %Find Number of entrys per event
            FE2 = FE(find(NE > 2));         %Filter events that are only single interactions
            NE2 = NE(find(NE > 2));         %Using NE2 for offsetting...
            
           for i = 1:NRead
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
    end

disp(['Found ',num2str(nfound),' events, Total events = ',num2str(ncoin)]);

% % Save Track Matrix, (optional)
% F = strfind(filename,'.dat');
% sname = [filename(1:F(1)-1),'_c.mat'];
% save(sname,'Res','chargeChange');
% disp(['Saved Seperated Events to, ' sname,' as a variable named Res{1,...}']);
    
    
    