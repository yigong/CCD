import numpy as np

def electron_hit(phi, azimuthal, geometry):
    '''
    Args:
        phi: electron scattered angle
        azimuthal: azimuthal angle electron. 0 is pointing the CCD
        geometry: a dict describes the geometry in mm.
    Return:
        True or false indicating whether the electron hits the CCD
    '''
    zWindow = geometry['verticle_distance']
    yEdge = geometry['distance_from_window']
    pxSize = 10.5/1000
    
    if 0.5*np.pi < azimuthal < 1.5*np.pi:
        xPixel = np.tan(azimuthal) * zWindow
        yPixel = np.sqrt(xPixel**2 + zWindow**2) * np.tan(phi)
        row = (yPixel - yEdge)/pxSize 
        col = xPixel/pxSize + 1754
        if 0 < row < 3508 and 0 < col < 3508:
            return (True, (int(row), int(col)))
        else:
            return (False, None)
        
    else:
        return (False, None)
    
    
    
#     x = geometry['verticle_distance']
#     
#     y = x * np.tan(azimuthal)
#     z = x /(np.cos(azimuthal) * np.tan(phi))
#     #print '(%s, %s) on yz plane.' % (y, z)
# 
#     
#     if z <= 0:
#         return False
#     
#     else:
#         CCD_center = (0, geometry['distance_from_window']+0.5*geometry['length'])
#         y_distance = y - CCD_center[0]
#         z_distance = z - CCD_center[1]
#         if abs(y_distance) < 0.5*geometry['width'] and abs(z_distance) < 0.5*geometry['length']:
# #             print '(%s, %s) on yz plane.' % (y, z)
# #             print 'y_distance = %.2f' % y_distance
# #             print 'z_distance = %.2f' % z_distance
#             return True
#         else:
#             return False
        
    
    # CCD_center is in (y,z)
            
        
    
    
    