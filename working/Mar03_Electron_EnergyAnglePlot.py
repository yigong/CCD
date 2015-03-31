from copy import deepcopy
Nitro_out = deepcopy(Nitro_10um)
Nitro_in = deepcopy(Nitro_10um)
keep = (Nitro_in.gammaEnergy > 1728)  & (Nitro_in.gammaEnergy < 1735)
Nitro_in.cutData(keep)
Nitro_out.cutData(~keep)
fig, ax = plt.subplots(figsize=(12,8))
angles = np.arccos(Nitro_in.cosPhi) * 180 / np.pi
ax.scatter(angles, Nitro_in.energyExit, color='b', lw=0, alpha=1, label='Within 10% bandwidth')
angles_org = np.arccos(Nitro_out.cosPhi) * 180 / np.pi
ax.scatter(angles_org, Nitro_out.energyExit, color='r', lw=0, alpha=0.1, label='Out of 10% bandwidth')

# angles_all = np.arccos(Nitro_10um.cosPhi) * 180 / np.pi
# ax.scatter(angles_all, Nitro_10um.energyExit, color='r', lw=0, alpha=0.5)

ax.set_xlabel('Electron scattering angles when exiting the scatter (deg)')
ax.set_ylabel('Electron energy (keV)')
ax.set_ylim((0, 1600))
ax.set_xlim((-2, 90))
#ax.legend()
fig.show()