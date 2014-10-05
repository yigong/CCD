from CCD.misc.unit_convertion import radius_to_degree
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, [ax0, ax1] = plt.subplots(1,2)

img0 = ax0.imshow(radius_to_degree(beta_min_map), origin='lower')
ax0.set_title('minimum beta in degree', fontsize='x-large')
ax0.get_xaxis().set_visible(False)
ax0.get_yaxis().set_visible(False)

divider = make_axes_locatable(ax0)
cax0 = divider.append_axes('right', size='5%', pad=0.1)
plt.colorbar(img0, cax=cax0)

img1 = ax1.imshow(energy_max_array, origin='lower')
ax1.set_title('maximum energy in keV', fontsize='x-large')
ax1.get_xaxis().set_visible(False)
ax1.get_yaxis().set_visible(False)

divider = make_axes_locatable(ax1)
cax1 = divider.append_axes('right', '5%', 0.1) 
fig.text(0.275,0.17,'Window')
fig.text(0.70,0.17,'Window')
fig.text(0.5,0.9,'1332.5 keV from Co-60', fontsize='x-large', ha='center')
plt.colorbar(img1, cax1)

figure()
hist(df.track_energy, bins=arange(0,1500,5), histtype='step')
xlabel('Energy (keV)', fontsize='x-large')
ylabel('counts', fontsize='x-large')
plt.show()



