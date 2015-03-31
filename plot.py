import matplotlib.pyplot as plt
import numpy as np
from CCD.io import load_pickle

def image2d(image_array, interactive_flag=True):
    from matplotlib.colors import Normalize
    
    cmap_Etrack = load_pickle('./cmap_Etrack.p') 
    norm_instance = Normalize(vmin=0, vmax=65535*5.546e-4)
    if interactive_flag == True:
        fig = plt.figure(figsize=(12,8))
        ax = fig.add_axes([0.1,0.1,0.6,0.8])
    else:
        fig, ax = plt.subplots()
    
    p = ax.imshow(image_array, cmap=cmap_Etrack, norm=norm_instance, interpolation='nearest', origin='lower')
    plt.show()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    return fig, ax
def format_coord(x, y):

    num_row = 100
    num_col = 100
    
    col = np.int(x + 0.5)
    row = np.int(y + 0.5)
   
    if col>=0 and col<=num_col and row>=0 and row<=num_row:
        z = 1
        
        return 'x=%1.4f, y=%1.4f, z=%1.4f' % (x,y,z)
    else:
        return 'x=%1.4f, y=%1.4f' % (x,y)
        