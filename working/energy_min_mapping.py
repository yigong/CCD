from CCD.io import read_file


threshold = 506 # the energy(keV) that gives an extrapolated range of 650 um.
ndim0 = ndim1 = 3506
energy_min = np.zeros([ndim0, ndim1])
rows, cols = np.where(energy_max_array < threshold)
energy_min[rows, cols] = energy_min_array[rows, cols]

row, cols = np.where(energy_min_array > threshold)
E_electron, dEdx = read_file('dEdx')
dxdE = 1e4/(dEdx*2329) # um/keV

