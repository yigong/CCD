import numpy as np

def degree_to_radius(input):
    return input/180.*np.pi

def radius_to_degree(input):
    return input/np.pi*180.