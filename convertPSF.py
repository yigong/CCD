from scipy.io import loadmat
psfTable = []
data = loadmat('psf.mat')
for i, p in enumerate(data['psf']):
    psfTable.append(p[30:91, 30:91])
    