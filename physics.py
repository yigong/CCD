def compton(gamma_energy, electron_angle):
    
    '''calculate for a compton event for given incoming gamma-ray energy and 
       scattered electron energy.
       
    Args:
        gamma_energy: incoming gamma-ray energy
        electron_angle: scattered electron angle in degree
           
    Return:
        result: a dict has following keys
            gamma_energy
            gamma_angle
            scattered_gamma_energy
            electron_energy
            electron_angle
            
    '''
    
    from sympy import sin, cos, cot, tan, solve, symbols
    import numpy as np
    
    electron_angle_radius = electron_angle*np.pi/180
    
    x = symbols('x')
    
    angle_tmp = solve( (1-cos(x))/511 + 1/gamma_energy - 
            (sin(x)/np.tan(electron_angle_radius)+cos(x))/gamma_energy, x)
    angle_tmp = np.array(angle_tmp)
    gamma_angle = angle_tmp[np.logical_and(angle_tmp>0, angle_tmp<np.pi)]
    
    gamma_angle = np.float(gamma_angle)
    
    scattered_gamma_energy = 1/((1-np.cos(gamma_angle))/511+1/gamma_energy)
    
    electron_energy = gamma_energy - scattered_gamma_energy
    
    result = {}
    result['gamma_energy'] = gamma_energy
    result['gamma_angle'] = gamma_angle
    result['scattered_gamma_energy'] = scattered_gamma_energy
    result['electron_energy'] = electron_energy 
    result['electron_angle'] = electron_angle

    return result
    