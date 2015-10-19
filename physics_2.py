import numpy as np

def Compton_electron_energy(gamma_energy, electron_angle, angle_type='deg'):
    mc2 = 511.
    gamma_energy = float(gamma_energy)
    if angle_type == 'deg':
        angles_rad = np.deg2rad(electron_angle)
    elif angle_type == 'rad':
        angles_rad = electron_angle
    part1 = ((mc2/gamma_energy + 1)/np.cos(angles_rad))**2 - 1
    return 2*mc2 / part1
