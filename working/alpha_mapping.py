from CCD.misc.unit_convertion import radius_to_degree


horizontal_distance = 4
lower_left = [horizontal_distance+0.5*10.5/1000, -1753.5*10.5/1000] # [dim0, dim1]
px_size = 10.5/1000

ndim0 = ndim1 = 3506
dim0_vector, dim1_vector = np.meshgrid(np.arange(ndim0), np.arange(ndim0), indexing='ij')
pos_window_coord = np.zeros([ndim0, ndim1])

pos_window_coord_dim0 = lower_left[0] + dim0_vector*px_size
pos_window_coord_dim1 = lower_left[1] + dim1_vector*px_size
# alpha map calculated from geometry
alpha_mean_map = np.arctan(-pos_window_coord_dim1/pos_window_coord_dim0)


beta_sample = np.linspace(10, 40, 7)*np.pi/180.
alpha_sigma_sample = np.array([9.80, 10.01, 10.52, 11.62, 12.96, 15.02, 18.15])

alpha_sigma_map = np.interp(beta_mean_map.flatten(), beta_sample, alpha_sigma_sample).reshape((ndim0, ndim1))
# for gaussian 95% --> 2 sigma each side
alpha_min_map = alpha_mean_map - 2*alpha_sigma_map
alpha_max_map = alpha_mean_map + 2*alpha_sigma_map
