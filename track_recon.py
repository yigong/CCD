import sys
from astropy.io import fits
from skimage.morphology import skeletonize
import numpy as np
from numpy import pi
from numpy.random import normal
from scipy.signal import convolve2d
from scipy.interpolate import RectBivariateSpline, UnivariateSpline
import pdb
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import re
import pickle
import traceback

def adjust_angles(angles_array):

    result =  angles_array - (abs(angles_array)>180) * np.sign(angles_array) * 360
    if np.any(result==180):
        result[result==180] = -180
    return result


def ridge_follow(fitsfile, outdir , plotflag=True, pickleflag=True):
    try:
    # if True:
        # initialization
        neighbr_matrix = np.ones([3,3])
        threshold = 0.2;
        num_step_back = 5;

        # find init
        file_idx = re.split(r'[/\.]', fitsfile)[-2]
        print 'recon #%s\n' %(file_idx)
        f = fits.open(fitsfile)
        image = f[0].data 
        true_pos = f[1].data
        image = image + normal(scale=0.025, size=np.shape(image))
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
        # for i, row_col in enumerate(zip(ends_row, ends_col)):
        #     print 'end %s: %s keV at (%s, %s)' %(i, ends_energy[i], ends_row[i], ends_col[i])
        # 'We pick end %s' %(end_idx)

        # step back 
        end_row_col= np.array((ends_row[end_idx], ends_col[end_idx])).astype('float')
        row_col = np.copy(end_row_col)
        thinned_copy= np.copy(thinned_image)
        for _ in xrange(num_step_back):
            thinned_copy[tuple(row_col)] = 0
            neighbors = thinned_copy[row_col[0]-1:row_col[0]+2, 
                                     row_col[1]-1:row_col[1]+2]
            row_col_rel = np.array(np.where(neighbors==1)) - 1
            row_col += row_col_rel[:,0]
        row_col_ridge = row_col.astype('float')
        rad2deg = lambda x : x*180./pi
        step_angle = rad2deg(np.arctan2(*(-row_col_rel)))
        alpha_guess = rad2deg(np.arctan2(*(end_row_col - row_col_ridge)))

        ## start ridge follow
        # initialization    
        num_transects = 60;
        row_num, col_num = np.shape(image)
        f_interp2d = RectBivariateSpline(range(row_num), range(col_num), image/16.) # 1/16 come from 0.25 grid size
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
        ridge_pos[0, :] = row_col_ridge
        ridge_angles = np.zeros(50)
        ridge_angles[0] = step_angle
        ridge_dEdx = np.zeros(50)
        ridge_centroid = np.zeros(50)
        transect_width = np.zeros(50)
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
            w = np.min(widths)
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

            transect_width[step_num-1] = w
            ridge_centroid[step_num] = centroid_value
            ridge_pos[step_num] = row_col_ridge
            ridge_angles[step_num] = step_angle
            ridge_dEdx[step_num-1] = dEdx
            threshold = .5
            if centroid_value < threshold/16.:
                break
        ridge_pos = ridge_pos[:step_num+1][:-10] # ignore last 10 ridge points
        ridge_angles = ridge_angles[:step_num+1]
        ridge_dEdx = ridge_dEdx[:step_num]
        ridge_centroid = ridge_centroid[1:step_num+1]
        transect_width = transect_width[:step_num]

        alpha_median = np.median(ridge_angles[1:-10])
        cov_matrix = np.cov(ridge_pos[1:, 1], ridge_pos[1:, 0])
        var_y = cov_matrix[1, 1]
        cov_xy = cov_matrix[0, 1]
        slope = var_y/cov_xy
        alpha_linearReg = adjust_angles(rad2deg(np.arctan2(-slope, -1)))

        # result dictionary
        result = dict()
        result['pos'] = ridge_pos
        result['angles'] = ridge_angles
        result['dEdx'] = ridge_dEdx
        result['centroid_value'] = ridge_centroid
        result['alpha_median'] = alpha_median
        result['alpha_linearReg'] = alpha_linearReg
        result['transect_width'] = transect_width

        # plot
        if plotflag:

            fig = Figure()
            canvas = FigureCanvasAgg(fig)
            ax = fig.add_subplot(111)
            ax.imshow(image, interpolation='nearest', origin='lower')
            ax.scatter(ridge_pos[:,1], ridge_pos[:,0], c='g', edgecolor='none')
            ax.plot(true_pos['col'], true_pos['row'], 'ro-', markersize=1.5, markeredgewidth=0)
            true_col_min = np.floor(min(true_pos['col']))
            true_col_max = np.ceil(max(true_pos['col']))
            true_row_min = np.floor(min(true_pos['row']))
            true_row_max = np.ceil(max(true_pos['row']))
            ax.set_xlim([true_col_min-2.5, true_col_max+2.5])
            ax.set_ylim([true_row_min-2.5, true_row_max+2.5])
            ax.set_title('threshold = %s || linear_reg: %.1f || median: %s  ' %(threshold, alpha_linearReg, alpha_median))
            canvas.print_figure('%s/%s_%s.png' %(outdir, 'ridge', file_idx))

            fig = Figure()
            canvas = FigureCanvasAgg(fig)
            ax = fig.add_subplot(111)
            ax.plot(ridge_dEdx, 'bs-', lw=2)
            ax.set_xlabel('ridge #')
            ax.set_ylabel('keV/um')
            ax.set_title('de/dx', fontsize=24)
            canvas.print_figure('%s/%s_%s.png' %(outdir, 'dEdx', file_idx))

            fig = Figure()
            canvas = FigureCanvasAgg(fig)
            ax = fig.add_subplot(111)
            ax.plot(np.arange(len(ridge_angles)), ridge_angles, 'bo-', lw=2)
            ax.set_xlabel('ridge #')
            ax.set_ylabel('alpha (deg)', color='blue')
            ax.set_ylim([-200, -120])
            for label in ax.get_yticklabels():
               label.set_color('blue')
            ax.set_title('alpha angles', fontsize=22)
            derivatives = np.diff(ridge_angles)
            ax1 = ax.twinx()
            ax1.plot(np.arange(len(derivatives))+0.5, derivatives, 'rs-', lw=1)
            ax1.set_ylabel('derivative of Alpha (deg)', color='red')
            ax1.set_ylim([-30, 80])
            for label in ax1.get_yticklabels():
                label.set_color('red')
            canvas.print_figure('%s/%s_%s.png' %(outdir, 'alpha', file_idx))

            fig = Figure()
            canvas = FigureCanvasAgg(fig)
            ax = fig.add_subplot(111)
            ax.plot(np.arange(len(transect_width)), transect_width, 'bo-', lw=2)
            ax.set_xlabel('ridge #')
            ax.set_ylabel('pix')
            ax.set_title('Cut width')
            canvas.print_figure('%s/%s_%s.png' %(outdir, 'width', file_idx))

        # output
        if pickleflag:
            pickle.dump(result, open('%s/%s.p' %(outdir, file_idx), 'wb'))
        return result

    except: 
   #  else:
        traceback.print_exc()
        tb = traceback.format_exc()
        open('%s/%s.err' %(outdir, file_idx), 'w').write(tb)
        return 0
