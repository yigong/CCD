num_tracks = length(recon_tracks);
reconTracks(num_tracks) = struct('track_position', [], 'track_original', [], 'track_1d', [],...
    'track_energy', [], 'track_thinned', [], 'ends_num', [], 'ends_pos', [], 'ends_idx', [],...
    'ends_energies', [], 'ridge_FWHMum', [], 'alpha', [], 'beta', []);

i = 1;
for j = 1: num_tracks
    if ~isempty(recon_tracks(j).alpha)
        reconTracks(i) = recon_tracks(j);
        i = i + 1;
    end
end
i = i - 1
reconTracks = reconTracks(1:i)