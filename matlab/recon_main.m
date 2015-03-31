clear
clc
%matlabpool local 4

load('/Users/Yigong/Google Drive/Research/CCD/data/segmented_tracks_doubleThickness')
%num_tracks = length(segmented_tracks_all);
num_tracks = 10000;
recon_tracks(num_tracks) = struct('track_position', [], 'track_original', [], 'track_1d', [],...
    'track_energy', [], 'track_thinned', [], 'ends_num', [], 'ends_pos', [], 'ends_idx', [],...
    'ends_energies', [], 'ridge_FWHMum', [], 'alpha', [], 'beta', []);

for i = 1:num_tracks
    
    track_init = segmented_tracks_all{i};
    recon_tmp = HybridTrack(track_init.track);
    if ~isfield(recon_tmp, 'err')
        recon_tracks(i).track_position = track_init.position;
        recon_tracks(i).track_original = track_init.track;
        recon_tracks(i).track_1d = track_init.track_1d;
        recon_tracks(i).track_energy = track_init.energy;

        recon_tracks(i).track_thinned = recon_tmp.EdgeSegments.thinnedImage;        
        recon_tracks(i).ends_num = recon_tmp.ends;
        recon_tracks(i).ends_pos = recon_tmp.EdgeSegments.coordinatesPix;
        recon_tracks(i).ends_idx = recon_tmp.EdgeSegments.chosenIndex;
        recon_tracks(i).ends_energies = recon_tmp.EdgeSegments.energiesKev;

        recon_tracks(i).ridge_FWHMum = recon_tmp.w;
        recon_tracks(i).alpha = recon_tmp.alpha;
        recon_tracks(i).beta = recon_tmp.beta;
        if mod(i, 100) == 0
            i
        end 
    end
    
end
%matlabpool close