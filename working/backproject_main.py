num_tracks = len(recon_tracks)
end_col_list = []
end_row_list = []

for i in range(num_tracks):
    alpha = recon_tracks[i]['alpha']
    end_idx = recon_tracks[i]['ends_idx'] - 1
    ends_num = recon_tracks[i]['ends_num']
    
    if ends_num != 1:
        end_pos = recon_tracks[i]['ends_pos'][end_idx]
        end_row = end_pos[0]
        end_col = end_pos[1]
        end_col_list.append(end_col)
        end_row_list.append(end_row)

    else:
        end_row = end_pos[0]
        end_row = end_pos[1]
        end_col_list.append(end_col)
        end_row_list.append(end_row)    
    EDGE = 4 # edge = 4 mm
    if 0<=alpha<=90 or 270<=alpha<=360:
        dim1 = (end_col - 1753 + 0.5) * 10.5/1000 + ((end_row+0.5)*10.5/1000 + EDGE) * np.tan(-1*alpha/180*pi)
    else:
        dim1 = None    
    recon_tracks[i]['backproject'] = dim1


    
    
  
    
    
    
    