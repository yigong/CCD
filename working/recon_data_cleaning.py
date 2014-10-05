from scipy.io import loadmat


EDGE = 4

tmp = loadmat(file_name)['recon_tracks']
recon_array = tmp[0]
dtype_names = recon_array.dtype.names
recon_tracks = []

for i in range(len(recon_array)): # loop the tracks
    recon_tmp = {} # the dict the save the track
    for j in range(len(dtype_names)): # loop the field
        key = dtype_names[j]
        recon_tmp[key] = np.squeeze(recon_array[i][key])

    lower_left = [np.min(recon_tmp['track_position'][0]), np.min(recon_tmp['track_position'][1])]
    recon_tmp['length'] = np.sum(recon_tmp['track_thinned'])
    recon_tmp['diffusivity'] = np.median(recon_tmp['ridge_FWHMum'])
    recon_tmp['ends_pos'] = recon_tmp['ends_pos'].astype(np.uint16) - 1 + lower_left# adjust for matlab 1 indexing
    recon_tmp['beta'] = recon_tmp['beta'].astype(np.float32)
    recon_tmp['ends_idx'] = recon_tmp['ends_idx'].astype(np.int16) - 1 # adjust for matlab 1 indexing
    recon_tmp['ends_num'] = recon_tmp['ends_num'].astype(np.int16)

    if recon_tmp['ends_num'] == 1:
        recon_tmp['ends_pos'] = [recon_tmp['ends_pos']]
        # to make sure recon_tmp['ends_pos']['ends_idx'] is a vector
    
    recon_tmp['px_num'] = np.sum(recon_tmp['track_1d'])
    
    end_idx = recon_tmp['ends_idx']
    # calculate back projected position on the window
    end_pos = recon_tmp['ends_pos'][end_idx]
    end_row = end_pos[0]
    end_col = end_pos[1]
    recon_tmp['init_segment'] = [end_row, end_col]
    alpha = recon_tmp['alpha']
    if 0<=alpha<=90 or 270<=alpha<=360:
        dim1 = (end_col - 1753 + 0.5) * 10.5/1000 + ((end_row+0.5)*10.5/1000 + EDGE) * np.tan(-1*alpha/180*pi)
    else:
        dim1 = None    
    recon_tmp['back-projection'] = dim1

    recon_tracks.append(recon_tmp)
    

    
    
    