from CCD.misc.unit_convertion import radius_to_degree
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, [ax0, ax1] = plt.subplots(1,2)
ndim0 = ndim1 = 3506
dim0_vector, dim1_vector = np.meshgrid(np.arange(ndim0), np.arange(ndim0), indexing='ij')
        
img0 = ax0.contour(dim1_vector,dim0_vector,  edata.alpha_mean_map, colors='k', origin='lower')
ax0.set_title('alpha calculated based on geometry', fontsize='x-large')
ax0.get_xaxis().set_visible(False)
ax0.get_yaxis().set_visible(False)

divider = make_axes_locatable(ax0)
cax0 = divider.append_axes('right', size='5%', pad=0.1)
plt.colorbar(img0, cax=cax0)

img1 = ax1.imshow(edata.alpha_max_map - edata.alpha_min_map, origin='lower')
ax1.set_title('alpha uncertainty (+/- give us 95% confidence)', fontsize='x-large')
ax1.get_xaxis().set_visible(False)
ax1.get_yaxis().set_visible(False)

divider = make_axes_locatable(ax1)
cax1 = divider.append_axes('right', '5%', 0.1) 
fig.text(0.275,0.17,'Window')
fig.text(0.70,0.17,'Window')
fig.text(0.5,0.9,'1332.5 keV from Co-60', fontsize='x-large', ha='center')
plt.colorbar(img1, cax1)
plt.show()




