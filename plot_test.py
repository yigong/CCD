from CCD.plot import format_coord



X = 10*np.random.rand(5,3)

figure();
imshow(X, cmap=plt.cm.jet, interpolation='nearest')
ax=gca()
ax.format_coord = format_coord
show()
 
