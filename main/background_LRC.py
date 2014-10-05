from tables import openFile

base_directory = '/Users/Yigong/Google Drive/Research/CCD/matlab/'

energy_top = np.array([])
energy_bottom = np.array([])


mat_name_tmp = base_directory + 'all_tracks.mat' 
file_tmp = openFile(mat_name_tmp)
energy_top = np.squeeze(file_tmp.root.energy_top[:])
energy_bottom = np.squeeze(file_tmp.root.energy_bottom[:])
tracks_top = np.squeeze(file_tmp.root.tracks_top[:])
tracks_bottom = np.squeeze(file_tmp.root.tracks_bottom[:])