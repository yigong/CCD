import numpy as np
import pandas as pd
from CCD.plot import image2d



df = pd.DataFrame(recon_tracks)

flag_length = np.logical_and(df.length >= length[0], df.length <= length[1])
if len(ends_num) == 1:
    flag_ends = df.ends_num == ends_num[0]
else:
    flag_ends = np.logical_and(df.ends_num >= ends_num[0], df.ends_num <= ends_num[1])
flag_diffusivity = np.logical_and(df.diffusivity >= diffusivity[0], df.diffusivity <= diffusivity[1])

flag_filter = np.logical_and(flag_length, flag_ends)
flag_filter = np.logical_and(flag_filter, flag_diffusivity)

df_filtered = df[flag_filter.values]

# define energy and alpha flag
flag_energy = pd.Series(np.zeros_like(flag_length), dtype=bool)
flag_energy.name = 'energy'

flag_alpha = pd.Series(np.zeros_like(flag_length), dtype=bool)
flag_alpha.name = 'alpha'

for i in df_filtered.index:
    alpha = recon_tracks[i]['alpha']
    energy = recon_tracks[i]['track_energy']
    init_idx = recon_tracks[i]['ends_idx']
    init_tmp = recon_tracks[i]['ends_pos'][init_idx]
    row_init = init_tmp[0]
    col_init = init_tmp[1]
    alpha_min = alpha_min_map[row_init, col_init]
    alpha_max = alpha_max_map[row_init, col_init]
    if alpha_min < alpha < alpha_max:
        flag_alpha[i] = True
     
    energy_min =  energy_min_map[row_init, col_init]
    energy_max =  energy_max_map[row_init, col_init]
    if energy_min < energy < energy_max:
        flag_energy[i] = True

flag_filter = np.logical_and(flag_filter, flag_alpha)
flag_filter = np.logical_and(flag_filter, flag_energy)

#flag_filter = pd.Series(np.ones_like(flag_length), dtype=bool) # include all tracks

df_filtered = df[flag_filter.values]


image_filtered = np.zeros((3506, 3506))
image_segmented_filtered = -1 * np.ones_like(image_filtered)
image_segmented_filtered = image_segmented_filtered.astype(np.int16)

for i in df_filtered.index:
    rows_in_image = recon_tracks[i]['track_position'][0]
    cols_in_image = recon_tracks[i]['track_position'][1]
    image_filtered[rows_in_image, cols_in_image] = recon_tracks[i]['track_1d'] + image_filtered[rows_in_image, cols_in_image]
#     row_min = recon_tracks[i]['track_position'][0]
#     col_min = recon_tracks[i]['track_position'][1]
#     track_original = recon_tracks[i]['track_original']
#     row_max = row_min + np.shape(track_original)[0]
#     col_max = col_min + np.shape(track_original)[1]
#     image_filtered[row_min:row_max, col_min:col_max] = track_original
    image_segmented_filtered[rows_in_image, cols_in_image] = int(i)





fig, ax = image2d(image_filtered)
text = fig.text(0.8,0.5,'track_info')

def on_mouse(event):
    event_row = int(event.ydata+0.5)
    event_col = int(event.xdata+0.5)
    
    segment_idx = (image_segmented_filtered[event_row, event_col])
    info = "Pixel row,col:{}, {}\n".format(event_row, event_col)
    if segment_idx != -1:
        track_info = recon_tracks[segment_idx]
        alpha_radius = -1 * track_info['alpha']/180.*np.pi
        end_idx = track_info['ends_idx'] 
        [end_row, end_col] = track_info['ends_pos'][end_idx]
        track_info['init_segment'] = [end_row, end_col] 
        print image_filtered[event_row, event_col]
        # draw an arrow at the end 
        arrow_length = 10
        ax.annotate("", xy=(end_col+arrow_length*np.sin(alpha_radius), end_row-arrow_length*np.cos(alpha_radius)),
                    xytext=(end_col, end_row), xycoords='data', arrowprops=dict(arrowstyle='->',color='g' ) )
        # add info
        
        if track_info['back-projection'] != None:
            print_tmp = """backproject: %(back-projection).1f mm
alpha: %(alpha).1f degree
energy: %(track_energy).1f keV
length: %(length)s pixels
diffusivity: %(diffusivity).1f um
init segment: %(init_segment)s
num of ends : %(ends_num)s
"""         % track_info
        else:
            print_tmp = """backproject: wrong direction
alpha: %(alpha).1f degree
energy: %(track_energy).1f keV
length: %(length)s pixels
diffusivity: %(diffusivity).1f um
init segment: %(init_segment)s
num of ends : %(ends_num)s
"""         % track_info
        
        info = info + print_tmp      
    else:
        #info = "Pixel row,col: %s, %s" % (event_row, event_col)
        print segment_idx
        
    text.set_text(info)
    fig.canvas.draw()

ax.set_title(file_name)
fig.show()
    
fig.canvas.mpl_connect('button_press_event', on_mouse)

