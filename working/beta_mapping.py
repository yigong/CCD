from CCD.physics import compton_electron_array

horizontal_distance = 4
lower_left = [horizontal_distance+0.5*10.5/1000, -1753.5*10.5/1000] # [dim0, dim1]
px_size = 10.5/1000

ndim0 = ndim1 = 3506
dim0_vector, dim1_vector = np.meshgrid(np.arange(ndim0), np.arange(ndim0), indexing='ij')
pos_window_coord = np.zeros([ndim0, ndim1])
pos_window_coord_dim0 = lower_left[0] + dim0_vector*px_size
pos_window_coord_dim1 = lower_left[1] + dim1_vector*px_size

pos_window_coord = np.sqrt(pos_window_coord_dim0**2 +\
                               pos_window_coord_dim1**2)


beta_min_map = (np.arctan(1/pos_window_coord)) 
beta_mean_map = (np.arctan(3.5/pos_window_coord))
beta_max_map = (np.arctan(6/pos_window_coord))

