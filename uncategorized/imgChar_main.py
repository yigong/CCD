# import the scripts
from imgChar_func import *

data_dir = '/Users/Yigong/Documents/beamCCD/data/'

# assign foldername, images want to be histogramed
foldername = data_dir + '04_24_14_160K/'
image_idx = range(23, 28)
row_idx = range(1766)
#
[baseline, noise] = histogram_img(folder=foldername, img_idx=image_idx, \
                                  plot_flag=True)
# f1, ax1 = plt.subplots()
#ax1.plot(baseline)


# 
# [baseline_row, noise_row] = img_char.histogram_row(folder=foldername, img_idx=image_idx, row_idx=row_idx, pflag=False)
# f1, ax1 = plt.subplots()
# ax1.plot(baseline_row)
# ax1.set_yscale('log')

plt.show()