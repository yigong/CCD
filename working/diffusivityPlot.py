import matplotlib.pyplot as plt
import numpy as np
import sys

BINS = np.arange(0, 40,01)



fig0, ax0 = plt.subplots()
ax0.hist(edataAfterPointBack.diffusivityMedian, bins=BINS, histtype='stepfilled',\
         normed='True', align='left',color='b', label='Co-60')
ax0.hold(True)
# beta gun data
#Sr90 = eData('../data/recon_tracks_beta_gun_8.mat')
(dummy1, dummy2, boxes) = ax0.hist(Sr90.diffusivityMedian, bins=BINS, histtype='step',\
                normed='True', align='left', color='r', label=r'beta emitter')
boxes[0].set_linewidth(2.0)
#boxes[0].set_linestyle('dashed')
ax0.legend()
ax0.set_xlim(0, 40)
ax0.set_ylim(0, 0.4)
ax0.set_xlabel('Diffusivity (um)')
ax0.set_ylabel('Normalized distribution')
ax0.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
fig0.show()
sys.exit('stop here')
# PB: point back
binsPB = np.arange(-20,20,0.5)
figPB, axPB = plt.subplots()
axPB.hist(edataAfterPointBack.back_projection, bins=binsPB, histtype='stepfilled',\
         align='left',color='g', label='Before width cut',alpha=0.5)
axPB.hold(True)
axPB.hist(edataAfterDiff.back_projection, bins=binsPB, histtype='stepfilled',\
         align='left',color='b', label='After width cut',alpha=1)
axPB.legend()
axPB.set_xlim(-22, 22)
axPB.set_ylim(0, 240)
axPB.set_xlabel('Pointing back position (mm)')
axPB.set_ylabel('Counts')

axPB.vlines([-2.75, 2.25], 0, 240, color='r',lw=1.5)
figPB.show()

# Ephi: phi vs. electron energy
# figEphi, axEphi = plt.subplots()
# axEphi.scatter(edataBeforeDiff.track_energy, edataBeforeDiff.phi, \
#                label='measurement', alpha=0.3, linewidths=0)
# axEphi.set_ylabel(r'$\hat{\phi}$ (degree)')
# axEphi.set_xlabel(r'$\hat{E_{e-}}$ (keV)')
# axEphi.plot(edataBeforeDiff.EeCompton[0], edataBeforeDiff.phiDegree, c='r', label='Compton Kinematics')
# axEphi.plot(edataBeforeDiff.EeCompton[1], edataBeforeDiff.phiDegree, c='r')
# axEphi.set_ylim(0, 90)
# axEphi.set_xlim(0, 1400)
# axEphi.legend()
# figEphi.show()





