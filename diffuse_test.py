from CCD.diffuse import diffuseAndPixelize
x = np.array([10, 20])
y = np.array([20, 40])
z = np.array([0, 100])
dE = np.array([10000, 10000])
psfTable = psf
pixelPlane = 'bottom'
image = diffuseAndPixelize(x, y, z, dE, psfTable, pixelPlane)
fig, ax = plt.subplots()
im = ax.imshow(image, origin='lower')
cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
fig.colorbar(im, cax=cax) 
plt.show()