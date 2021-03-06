import sys
from astropy.io import fits
from skimage.morphology import skeletonize
import numpy as np
from numpy import pi
from numpy.random import normal
from scipy.signal import convolve2d
from scipy.interpolate import RectBivariateSpline, UnivariateSpline
import pdb
from pylab import *
import re

def adjust_angles(angles_array):

    result =  angles_array - (abs(angles_array)>180) * np.sign(angles_array) * 360
    if any(result==180):
        result[result==180] = -180
    return result


def find_init(fits_file, plot_flag=True):
#    args = sys.argv
#    fits_file = args[1] 
    
    # initialization
    neighbr_matrix = np.ones([3,3])
    threshold = 0.2;
    num_step_back = 5;

    # find init
    file_idx = re.split(r'[/\.]', fits_file)[-2]
    f = fits.open(fits_file)
    image = f[0].data 
    true_pos = f[1].data
    image = image + normal(scale=0.03, size=np.shape(image))
    binary_image = (image > threshold).astype(np.int8)    # binary image
    thinned_image = skeletonize(binary_image).astype(np.int8)
    neighbr_num_image = convolve2d(thinned_image, neighbr_matrix, 'same')
    neighbr_num_thinned_image = neighbr_num_image * thinned_image
    ends_row, ends_col = np.where(neighbr_num_thinned_image == 2) 
    sum_matrix = np.ones([6,6])
    sum_matrix[0, (0,1,-2,-1)] = 0
    sum_matrix[1, (0,-1)] = 0
    sum_matrix[-2, (0,-1)] = 0
    sum_matrix[-1, (0,1,-2,-1)] = 0
    ends_energy = convolve2d(image, sum_matrix, 'same')[ends_row, ends_col]
    end_idx = np.argmin(ends_energy)
    for i, row_col in enumerate(zip(ends_row, ends_col)):
        print 'end %s: %s keV at (%s, %s)' %(i, ends_energy[i], ends_row[i], ends_col[i])
    'We pick end %s' %(end_idx)

    # step back 
    end_row_col= np.array((ends_row[end_idx], ends_col[end_idx]))
    row_col = np.copy(end_row_col)
    thinned_copy= np.copy(thinned_image)
    for _ in xrange(num_step_back):
        thinned_copy[tuple(row_col)] = 0
        neighbors = thinned_copy[row_col[0]-1:row_col[0]+2, 
                                 row_col[1]-1:row_col[1]+2]
        row_col_rel = np.array(np.where(neighbors==1)) - 1
        row_col += row_col_rel[:,0]
    RF_row_col = row_col
    rad2deg = lambda x : x*180./pi
    stepAngle = rad2deg(np.arctan2(*(-row_col_rel)))
    alphaGuess = rad2deg(np.arctan2(*(end_row_col - RF_row_col)))

    # plot
    # if plot_flag:
       #  figure()
       #  imshow(binary_image, interpolation='nearest', origin='lower')
       #  figure()
       #  imshow(thinned_image, interpolation='nearest', origin='lower')
       #  scatter(*(end_row_col[::-1]))
       #  scatter(*(RF_row_col[::-1]))
       #  show()
    # assign values to output dict
    result = dict() 
    result['end_row_col'] = end_row_col
    result['RF_row_col'] = row_col
    result['step_angle'] = stepAngle
    result['alpha_guess'] = alphaGuess
    result['image'] = image
    result['true_pos'] = true_pos
    result['file_idx'] = file_idx
    return result 

def ridge_follow(init, plot_flag=True, plot_save_dir=None):

    # initialization    
    num_transects = 60;
    image = init['image']
    row_num, col_num = np.shape(image)
    alpha_guess = init['alpha_guess']
    step_angle = init['step_angle']
    f_interp2d = RectBivariateSpline(range(row_num), range(col_num), image/16.) # 1/16 come from 0.25 grid size
    row_col_ridge = init['RF_row_col'].astype('float')
    true_pos = init['true_pos'] 
    # pre-compute cuts
    col_0deg = np.arange(-5, 5.001, 0.25)
    row_0deg = np.zeros_like(col_0deg)
    deg2rad = lambda x : x/180.*pi
    angles_rad = deg2rad(np.arange(180))
    cos_v, length_v = np.meshgrid(np.cos(angles_rad), col_0deg, indexing='ij')
    sin_v, _        = np.meshgrid(np.sin(angles_rad), col_0deg, indexing='ij')
    row_rotations = length_v * sin_v
    col_rotations = length_v * cos_v

    # ridge follow
    ridge_pos = np.zeros([50, 2])
    ridge_pos[0, :] = init['RF_row_col']
    ridge_angles = np.zeros(50)
    ridge_angles[0] = step_angle
    ridge_dEdx = np.zeros(50)
    ridge_centroid = np.zeros(50)
    for step_num in xrange(1, 50):
        transects_angle_center = adjust_angles(step_angle+90) 
        transects_angle_array = adjust_angles(np.arange(int(transects_angle_center-num_transects/2), int(transects_angle_center+num_transects/2)+0.1)).astype(int)
        row_transects = row_rotations[transects_angle_array] + row_col_ridge[0]
        col_transects = col_rotations[transects_angle_array] + row_col_ridge[1]
        transects_shape = np.shape(col_transects)
        energy_transects = f_interp2d(row_transects.flatten(), col_transects.flatten(), grid=False).reshape(transects_shape)
        widths = np.ones_like(transects_angle_array) * 100
        for i, (angle, dE) in enumerate(zip(transects_angle_array, energy_transects)):
            spline = UnivariateSpline(col_0deg, dE - max(dE)/2.)
            roots = spline.roots()
            if len(roots) == 2:
                widths[i] = roots[1] - roots[0]
        # useful arrays: transects_angle_array, centroid_indices, centroid_values, width
        transect_chosen_index = np.argmin(widths)
        transect_chosen = energy_transects[transect_chosen_index]
        transect_angle_chosen = transects_angle_array[transect_chosen_index]
        centroid_value = np.max(transect_chosen)
        centroid_index = np.argmax(transect_chosen)
        centroid_pos   = [row_transects[transect_chosen_index, centroid_index], col_transects[transect_chosen_index, centroid_index]]
        step_angle = (transect_angle_chosen-90)
        step_vector = 0.25*np.array([np.sin(deg2rad(step_angle)), np.cos(deg2rad(step_angle))])
        row_col_ridge = centroid_pos + step_vector
        # row_col_ridge = row_col_ridge + step_vector
        dEdx = np.sum(transect_chosen) / (0.25 * 10.5) # in keV/um

        ridge_centroid[step_num] = centroid_value
        ridge_pos[step_num] = row_col_ridge
        ridge_angles[step_num] = step_angle
        ridge_dEdx[step_num-1] = dEdx
        threshold = .5
        if centroid_value < threshold/16.:
            break
    ridge_pos = ridge_pos[:step_num+1]
    ridge_angles = ridge_angles[:step_num+1]
    ridge_dEdx = ridge_dEdx[:step_num]
    ridge_centroid = ridge_centroid[1:step_num+1]

    # plot
    if plot_flag:
        file_idx = init['file_idx']
        figure()
        imshow(image, interpolation='nearest', origin='lower')
        scatter(ridge_pos[:,1], ridge_pos[:,0], c='g')
        plot(true_pos['col'], true_pos['row'], 'ro-', markersize=1.5, markeredgewidth=0)
        xlim([-0.5, col_num-0.5])
        ylim([-0.5, row_num-0.5])
        title('threshold = %s' %(threshold))
        savefig('%s/%s_%s.png' %(plot_save_dir, 'ridge', file_idx))

        figure()
        plot(ridge_dEdx, 'bs-', lw=2)
        xlabel('ridge #')
        title('dE/dx', fontsize=24)
        savefig('%s/%s_%s.png' %(plot_save_dir, 'dEdx', file_idx))

        figure()
        plot(np.arange(len(ridge_angles)), ridge_angles, 'bo-', lw=2)
        xlabel('ridge #')
        ylabel('Alpha (deg)', color='blue')
        ylim(ymin = min(ridge_dEdx)-30)
        title('Alpha angles', fontsize=22)
        derivatives = np.diff(ridge_angles)
        ax1 = twinx()
        ax1.plot(np.arange(len(derivatives))+0.5, derivatives, 'rs-', lw=1)
        ylabel('derivative of Alpha (deg)')
        ylim(ymax = max(derivatives)+30)
        savefig('%s/%s_%s.png' %(plot_save_dir, 'alpha', file_idx))
        # show()

    alpha_median = median(ridge_angles[1:])
    cov_matrix = np.cov(ridge_pos[1:, 1], ridge_pos[1:, 0])
    var_y = cov_matrix[1, 1]
    cov_xy = cov_matrix[0, 1]
    slope = var_y/cov_xy
    alpha_linearReg = adjust_angles(rad2deg(np.arctan2(slope, 1))) - 180

    # result dictionary
    result = dict()
    result['pos'] = ridge_pos
    result['angles'] = ridge_angles
    result['dEdx'] = ridge_dEdx
    result['centroid_value'] = ridge_centroid
    result['alpha_median'] = median(ridge_angles)
    result['alpha_linearReg'] = alpha_linearReg

    1/0
    return result
    



if __name__ == '__main__':
    init = find_init()
    ridge = ridge_follow(init)
