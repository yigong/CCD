def convert():
    from glob import glob
    import os
    
    filenames = glob('/Users/Yigong/Documents/beamCCD/data/07_25_14_Y88/Y88*.fit')
    for filename in filenames:
        suffix = filename.split('image')[-1]
        os.rename(filename, 'image_%s' % suffix)
    
    return