from PyBeamDiag.physics_2 import Compton_electron_energy
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import pickle
from matplotlib.colors import Normalize

row_1d = np.arange(3700)
col_1d = np.arange(3700)
row_vec, col_vec = np.meshgrid(row_1d, col_1d, indexing='ij')
y = 2000 + 10.5/2 + row_vec * 10.5
x = -3700./2*10.5 + 10.5/2 + col_vec * 10.5
z = np.zeros_like(x)
z0 = 3500. # in um

mc2 = 511.
cosPhi = y/(x**2 + y**2 + (z-z0)**2)**(0.5)
Phi = np.arccos(cosPhi)
Ee = Compton_electron_energy(1732., Phi, angle_type='rad')
tanA = np.abs(x)/y
tanB = (z0-z)/(x**2 + y**2)**(0.5)
dz = 1000. # thickness of the strip in um
tanB_max = (z0 + 0.5*dz - z)/(x**2 + y**2)**(0.5)
tanB_min = (z0 - 0.5*dz - z) / (x**2 + y**2)**(0.5)
beta = np.deg2rad(np.arange(5, 46, 5))
dAlpha = np.deg2rad(np.array([13.1, 13.3, 14.1, 15.9, 17.8, 20.8, 25.4, 32.2, 38.0]))
dA = np.interp(np.arctan(tanB), beta, dAlpha)
dB = np.arctan(tanB_max) - np.arctan(tanB_min)

c1 = cosPhi * (1 + 2*mc2/Ee)**(0.5)
Er = mc2 / (c1 - 1)
fig = Figure()
canvas = FigureCanvasAgg(fig)
ax = fig.add_subplot(111)
im = ax.imshow(Er, interpolation='nearest', origin='lower', norm=Normalize())
fig.colorbar(im)
ax.set_xlabel('x (um)')
ax.set_ylabel('y (um)')
ax.set_title('Er')
canvas.print_figure('Er.png')

A1 = mc2 * c1 / (c1 - 1)**2
Er_FWHM = A1 * (tanA**2 * dA**2 + tanB**2 * dB**2)**(0.5)
fig = Figure()
canvas = FigureCanvasAgg(fig)
ax = fig.add_subplot(111)
im = ax.imshow(Er_FWHM, interpolation='nearest', origin='lower', norm=Normalize(vmin=10, vmax=100))
fig.colorbar(im)
ax.set_xlabel('x (um)')
ax.set_ylabel('y (um)')
ax.set_title('FWHM Er')
canvas.print_figure('Er_FWHM.png')
result = dict()
result['Er'] = Er
result['Er_FWHM'] = Er_FWHM
# pickle.dump(result, open('Er.p', 'wb'))


