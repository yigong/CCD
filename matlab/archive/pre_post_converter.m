GAIN = [5.609e-4, 5.4834e-4];
energy = [energy_bottom*GAIN(1), energy_top*GAIN(2)];

num_tracks_bottom = length(tracks_bottom);
num_tracks_top = length(tracks_top);
num_tracks = num_tracks_bottom + num_tracks_top;

tracks_tmp = [tracks_bottom{:}, tracks_top{:}];
tracks(num_tracks) = struct;

gain = GAIN(1);
for i = 1:num_tracks
    
    if i == num_tracks_bottom+1
        gain = GAIN(2);
        'now processing top'
    end
    
    tracks(i).img = tracks_tmp(i).img * gain;
    tracks(i).E = tracks_tmp(i).E * gain;
    tracks(i).x = tracks_tmp(i).x;
    tracks(i).y = tracks_tmp(i).y;


end
