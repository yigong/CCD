import numpy as np

def electron_hit(theta, phi, geometry):
    '''
    Args:
        theta: electron scattered angle
        phi: azimuthal angle electron. 0 is pointing the CCD
        geometry: a dict describes the geometry in mm.
    Return:
        True or false indicating whether the electron hits the CCD
    '''
    
    x = geometry['verticle_distance']
    
    y = x * np.tan(phi)
    z = x /(np.cos(phi) * np.tan(theta))
    #print '(%s, %s) on yz plane.' % (y, z)

    
    if z <= 0:
        return False
    
    else:
        CCD_center = (0, geometry['distance_from_window']+0.5*geometry['length'])
        y_distance = y - CCD_center[0]
        z_distance = z - CCD_center[1]
        if abs(y_distance) < 0.5*geometry['width'] and abs(z_distance) < 0.5*geometry['length']:
#             print '(%s, %s) on yz plane.' % (y, z)
#             print 'y_distance = %.2f' % y_distance
#             print 'z_distance = %.2f' % z_distance
            return True
        else:
            return False
        
    
    # CCD_center is in (y,z)
            
        
    
    
    