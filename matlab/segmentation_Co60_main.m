warning off
base_dir = './'
base_dir = './'
flist = dir([base_dir, 'image_5_corrected.fit']);

num_tmp = size(flist);
num_file = num_tmp(1)
dirlist{num_file, 1} = [];

tracks = {};
energy = [];

tracks_tmp = {};
energy_tmp = [];
THRESHOLD = [ 0.125 ,0.181 ];
GAIN = [5.6090e-4, 5.4834e-4];
j = 0;

for i = 1: num_file
 
    img = fitsread(fullfile(base_dir, dirlist{i}, flist(i).name));
    
    [T, E] = CCDsegment4(img, 6, 'Mode', 'postcal', 'Gain', GAIN);
    
    sprintf('%s tracks are segmented from  %sth image', int2str(length(E)), int2str(i))

    tracks_tmp = [tracks_tmp, T];
    energy_tmp = [energy_tmp, E];
%     
%     if  length(energy_tmp) > 10000
%         j = j + 1;
%         tracks = [tracks, tracks_tmp];
%         energy = [energy, energy_tmp];
%         sprintf('segmenting image %s out of %s.', int2str(i), int2str(num_file) )
%         sprintf('we have %s tracks now.', int2str(length(energy)))
%         
%         save_name = sprintf('./tracks_%s.mat', int2str(j))
%         save(save_name, 'tracks_tmp', 'energy_tmp', '-v7.3')
%         
%         tracks_tmp = {};
%         energy_tmp = [];
%         
%     
%     end
%             tracks = [tracks, tracks_tmp];
    energy = [energy, energy_tmp];
    tracks = [tracks, tracks_tmp];

end

save('./tracks.mat', 'energy', 'tracks', '-v7.3')
