from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
import matplotlib.cm as cmx

jet = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=1200)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

Emin = edata.energy_min_map
Emax = edata.energy_max_map
figE, [axEmin, axEmax] = plt.subplots(1,2)
imgEmin = axEmin.imshow(Emin, origin='lower', norm=cNorm)
imgEmax = axEmax.imshow(Emax, origin='lower', norm=cNorm)

axEmin.get_xaxis().set_visible(False)
axEmin.get_yaxis().set_visible(False)
axEmax.get_xaxis().set_visible(False)
axEmax.get_yaxis().set_visible(False)
# divider = make_axes_locatable(axEmax)
# cax0 = divider.append_axes('right', size='5%', pad=0.1)
# plt.colorbar(imgEmax, cax=cax0)
# 
# divider = make_axes_locatable(axEmin)
# cax1 = divider.append_axes('right', size='5%', pad=0.1)
# plt.colorbar(imgEmin, cax=cax1)




ndim0 = 3506
LV = [10, 100, 300, 500, 700, 900, 1100]
row, col = np.meshgrid(np.arange(ndim0), np.arange(ndim0), indexing='ij')
CS_Emin = axEmin.contour(col, row, Emin, colors='k', linewidths=2, levels=LV, origin='lower')
CS_Emax = axEmax.contour(col, row, Emax, colors='k', linewidths=2, levels=LV, origin='lower')
Emin_locations = [(300, 600), (450, 1200), (600, 1800), (800, 2700), (1500, 2900)]
Emax_locations = [(200, 350), (350, 750), (500, 1050), (650, 1400), (800, 1950), (1500, 1500)]
plt.clabel(CS_Emin, fontsize=20, inline=1, fmt='%.d', manual=Emin_locations, inline_spacing=1)
plt.clabel(CS_Emax, fontsize=20, inline=1, fmt='%.d', manual=Emax_locations, inline_spacing=1)
figE.show()