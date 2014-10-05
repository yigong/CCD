from CCD.physics import compton_electron_array
from CCD.misc.unit_convertion import radius_to_degree
from mpl_toolkits.axes_grid1 import make_axes_locatable

horizontal_distance = 4
lower_left = [horizontal_distance+0.5*10.5/1000, -1753.5*10.5/1000] # [dim0, dim1]
px_size = 10.5/1000

ndim0 = ndim1 = 3506
dim0_vector, dim1_vector = np.meshgrid(np.arange(ndim0), np.arange(ndim0), indexing='ij')
pos_window_coord = np.zeros([ndim0, ndim1])
pos_window_coord_dim0 = lower_left[0] + dim0_vector*px_size
pos_window_coord_dim1 = lower_left[1] + dim1_vector*px_size

# pos_window_coord = np.sqrt(pos_window_coord_dim0**2 +\
#                                pos_window_coord_dim1**2)
opposite = np.sqrt(3.5**2 + pos_window_coord_dim1**2)
adjacent = pos_window_coord_dim0 
phi_min_map = (np.arctan((opposite-2.5)/adjacent)) 
phi_mean_map = (np.arctan(opposite/adjacent))
phi_max_map = (np.arctan((opposite+2.5)/adjacent))

energy_min_map = compton_electron_array(1173.2, phi_max_map.flatten()).reshape((ndim0, ndim1))
energy_mean_map = compton_electron_array(1332.5, phi_mean_map.flatten()).reshape((ndim0, ndim1))
energy_max_map = compton_electron_array(1332.5, phi_min_map.flatten()).reshape((ndim0, ndim1))
