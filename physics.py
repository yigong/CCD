import numpy as np

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

def compton_gamma(gamma_energy, gamma_scattered_angle):
    electron_rest_mass = 511
    
    tmp = (1 - np.cos(gamma_scattered_angle))/electron_rest_mass
    
    gamma_scattered_energy = 1/(tmp+1/gamma_energy)
    
    ratio = gamma_energy/gamma_scattered_energy
    
    cot_electron_angle = (ratio - np.cos(gamma_scattered_angle))/np.sin(gamma_scattered_angle)
    electron_scattered_angle = np.arctan(1/cot_electron_angle)
    
    electron_energy = gamma_energy-gamma_scattered_energy
    
    
    
    result = {}
    result['gamma_energy'] = gamma_energy
    result['gamma_scattered_angle'] = gamma_scattered_angle
    result['gamma_scattered_energy'] = gamma_scattered_energy
    result['electron_energy'] = electron_energy
    result['electron_scattered_angle'] = electron_scattered_angle
    
    return result
    
    

def Klein_Nishina_pdf(gamma_energy,num_samples):
    ''' sampling Klein Nishina from cos(theta)
    
    Args:
        number of samples
    Return:
        an array of sampled pdf
        
    '''
    electron_rest_mass = 511.
    alpha = gamma_energy/electron_rest_mass
    
    cos_polar = np.linspace(1,-1, num_samples )
    r = 1/(1+alpha*(1-cos_polar))
    # r is the ratio between scattered photon energy and initial energy 
    pdf_tmp = r**2*(r + 1/r - 1 + cos_polar**2)/2
    pdf = pdf_tmp/np.sum(pdf_tmp)
    
    return pdf, cos_polar
    
    
    