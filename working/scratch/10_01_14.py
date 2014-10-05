figure();
hist(edata.gammaEnergy, bins=linspace(0, 1400, 1401))
xlabel('Gamma-ray energy (keV)', fontsize='xx-large')
ylabel('Counts', fontsize='xx-large')
vlines([1332, 1173], 0, 2, color='r', alpha=0.5, linestyles='dashed')


figure();
scatter(edata.track_energy, np.array(edata.phi)*180/np.pi)
xlabel('Electron energy (keV)',  fontsize='xx-large')
ylabel('$\phi$(degree)', fontsize='xx-large')


show()