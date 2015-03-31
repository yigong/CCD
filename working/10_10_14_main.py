from CCD.working.track_analysis import eData
from CCD.working.PositionSampling import PositionSampling


figure();
hist(edata.gammaEnergy, bins=linspace(0, 1400, 141))
xlabel('Gamma-ray energy (keV)', fontsize='xx-large')
ylabel('Counts', fontsize='xx-large')
vlines([1332, 1173], 0, 4, color='r', alpha=0.5, linestyles='dashed')



figure();
scatter(edata.track_energy, np.array(edata.phi)*180/np.pi)
xlabel('Electron energy (keV)',  fontsize='xx-large')
ylabel('$\phi$(degree)', fontsize='xx-large')


show()




fig, ax = plt.subplots()
edata_30s = eData(filename='./recon_tracks_73min.mat')
ax.hist(edata_30s.ends_num, bins=range(1,10), normed=True, histtype='step', \
        align='left',label='30s exposure')

edata_60s = eData(filename='./recon_tracks_45000.mat')
ax.hist(edata_60s.ends_num, bins=range(1,10), normed=True, histtype='step', \
        align='left',label='60s exposure', color='red')

ax.set_xlabel('Number of ends', fontsize='xx-large')
ax.set_ylabel('Normalized Distribution', fontsize='xx-large')
legend()
plt.show()