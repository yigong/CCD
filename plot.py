import matplotlib.pyplot as plt
import numpy as np

def image2d(image_array, color_map=plt.cm.gray):
    fig, ax = plt.subplots()
    ax.imshow(image_array, cmap=color_map, interpolation='nearest')
    plt.show()