import sys
from glob import glob
from CCD.track_recon import find_init, ridge_follow

def main():
    args = sys.argv
    f_dir = args[1] # fits directory
    p_dir = args[2] # plot directory

    f_files = glob('%s/*.fits' %(f_dir))
    for f in f_files:
        init = find_init(f)
        ridge = ridge_follow(init, plot_save_dir=p_dir)

if __name__ == '__main__':
    main()
