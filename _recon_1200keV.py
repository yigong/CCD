from glob import glob
import sys
from PyBeamDiag.track_recon import ridge_follow
import multiprocessing
from functools import partial

def main():
    args = sys.argv
    fits_dir = args[1]
    pickle_dir = args[2]
    f_files = glob('%s/*.fits' %(fits_dir))
    for f in f_files:
        ridge_follow(f, pickle_dir)
#     partial_func = partial(ridge_follow, out_dir=pickle_dir)
#     pool = multiprocessing.Pool(15)
#     pool.map(partial_func, f_files)

if __name__ == '__main__':
    main()
