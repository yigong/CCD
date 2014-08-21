import numpy as np

def extract_tracks(image_stack, gain, threshold):
    '''extract tracks from the CCD image
    
    Args:
        image_stack: median subtracted image stack
        gain: ~
        threshold: threshold for the segmentation, in ADC value
    
    Returns:
        a list of tracks and its position
    
    
    '''
    
    from skimage.morphology import binary_dilation
    from skimage.measure import label
    
    neighbour_array = np.ones((3,3))
    track_list = []
    length_list = []
    energy_list = []
    position_list = []
    for image_idx in range(len(image_stack)):
        image_temp = image_stack[image_idx, ...]
        
        image_binary = image_temp > threshold
        image_dilated = binary_dilation(image_binary, neighbour_array)
        #image_dilated = image_binary
        image_segmented, num_label = label(image_dilated, neighbors=8, 
                                           return_num=True)
        print '%s tracks in image # %s (out of %s)' % (num_label, image_idx+1, len(image_stack))
        
        for track_idx in range(1, num_label):
            if track_idx%1000 == 0:
                print track_idx
            rows, cols = np.where(image_segmented == track_idx)
            row_min = np.min(rows)
            row_max = np.max(rows)
            col_min = np.min(cols)
            col_max = np.max(cols)
            
            track_temp = image_temp[row_min: row_max+1, col_min: col_max+1]
            #track_list.append(track_temp * gain)
            
            length_list.append(len(rows))
            
            pixels_lit_up = image_temp[rows, cols]
            energy_list.append(np.sum(pixels_lit_up) * gain)
            
            #position_list.append((row_min, col_min))
            
            #1/0
            
#             single_dict = {'track' : track_array*gain,
#                           'energy' : np.sum(track_array)*gain,
#                           'position' : (row_min, col_min)
#                           }
#             # the position is the (row, col) of the lower left pixel in the track.
#             track_dict[str(image_idx)] = single_dict
            
    return (track_list, energy_list, length_list, position_list)