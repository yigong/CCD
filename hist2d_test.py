from numpy.random import rand
a = rand(100)
b = 9 + rand(100)
H, c, d = np.histogram2d(a, b, bins=np.arange(11))
fig, ax = plt.subplots()
imgPlot = ax.imshow(H, origin='lower', interpolation='nearest')
plt.colorbar(imgPlot)
plt.show()