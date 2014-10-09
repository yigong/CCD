fig, ax = plt.subplots()
ax.hist(edata.ends_num, bins=range(1,10), normed=True, histtype='step', align='left')
ax.set_xlabel('Number of ends', fontsize='xx-large')
ax.set_ylabel('Normalized Distribution', fontsize='xx-large')
plt.show()