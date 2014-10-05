from IPython import get_ipython
import pandas as pd

file_name = '/Users/Yigong/Google Drive/Research/CCD/working/recon_tracks_45000.mat'

length = [0, np.inf] 
ends_num =  [2]
diffusivity = [0, 15]


ipython = get_ipython()
scripts_to_run = ['beta_mapping.py',
                  'alpha_mapping.py',
                  'phi_mapping.py',
                  'recon_data_cleaning.py']

for script in scripts_to_run:
    print 'running %s' % script
    ipython.magic(r'%%run -i %s' % script)
df = pd.DataFrame(recon_tracks)