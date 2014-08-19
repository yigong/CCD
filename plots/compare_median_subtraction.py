import matplotlib.pyplot as plt
from CCD.io import load_pickle

energy_list = load_pickle('./compare_median_subtraction.p_plot')
fig, (ax0, ax1) = plt.subplots(ncols=2, sharey=True, figsize=(14,6))


ax0.hist(energy_list['pixel-by-pixel'], bins=np.arange(4000, 200000, 100))
ax0.set_xlabel('track in ADC')
ax0.set_ylabel('counts')
ax0.set_ylim(0,300)
ax0.set_title('pixel-by-pixel median subtraction algorithm')

ax1.hist(energy_list['self'], bins=np.arange(4000, 200000, 100))
ax1.set_xlabel('track in ADC')
ax1.set_title('self median subtraction algorithm')

plt.show()


