flag = np.logical_and(edataOriginal.track_energy>600, edataOriginal.track_energy<620)
edataOriginal.cutData(flag)
edataOriginal.plotInteractive()