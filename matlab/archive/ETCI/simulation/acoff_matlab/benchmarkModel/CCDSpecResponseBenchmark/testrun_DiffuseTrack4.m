dE = []; 
progressbar(0);

load ../../../../psft_650um.mat

for i=1:length(Event); 
%     out{i} = Geant4TrackHandling(Event{i}.trackM,'coinc','nonoise','zcentered','psft',psft); 
    out{i} = Geant4TrackHandling(Event{i}.trackM,'coinc','psft',psft);
    if ~isempty(out{i}) && ~isempty(out{i}.E); 
        dE(i) = out{i}.E(1) - Event{i}.Edep(1); 
    else
        dE(i) = nan; 
    end
    progressbar(i/length(Event));
end