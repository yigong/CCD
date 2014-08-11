import numpy as np
import matplotlib.pyplot as plt

def myplot(*args, **kwargs):
    f1, ax1 = plt.subplots()
    ax1.plot(*args, **kwargs)
    plt.show()
    return f1, ax1
    
def myhist(*args, **kwargs):
    f1, ax1 = plt.subplots()
    ax1.hist(*args, **kwargs)
    plt.show()
    return f1, ax1