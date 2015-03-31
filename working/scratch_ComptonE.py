import numpy as np
import matplotlib.pyplot as plt

# plot electron energy vs. electron scattering angle
Eg = 1332.
mc2 = 511. # electron rest mass energy
#cosPhiMin = np.sqrt(Eg*(Eg+2*mc2)) # find the Maximum scattering angle

phiDegree = np.linspace(0, 90, 100)
phi = phiDegree/180.*np.pi
Ee = 2*mc2/(((mc2/Eg+1)/np.cos(phi))**2 - 1)
fig0, ax0 = plt.subplots()
ax0.plot(phiDegree, Ee)
ax0.set_xlabel(r'$\hat{\phi}$', fontsize='xx-large')
ax0.set_ylabel('$\hat{E_{e-}}$', fontsize='xx-large')
            
plt.show()

