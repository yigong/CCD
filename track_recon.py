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
        if isinstance(result, (float, int)):
            result = -180
        else:
            result[result==180] = -180
    return result


def ridge_follow(fitsfile, outdir , plotflag=True, pickleflag=True):
    try:
    # if True:
        # initialization
        neighbr_matrix = np.ones([3,3])
        binThreshold = 0.5;
        num_step_back = 4;

        # find init
        file_idx = re.split(r'[/\.]', fitsfile)[-2]
        print 'recon #%s\n' %(file_idx)
        f = fits.open(fitsfile)
        image = f[0].data 
        true_pos = f[1].data
        image = image + normal(scale=0.025, size=np.shape(image))
        binary_image = (image > binThreshold).astype(np.int8)    # binary image
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
            if len(row_col_rel) > 1:
                row_col_rel = row_col_rel[:,0]
            row_col += row_col_rel
        ridge_row_col = row_col.astype('float')
        rad2deg = lambda x : x*180./pi
        stepD = rad2deg(np.arctan2(*(-row_col_rel)))
        alpha_guess = rad2deg(np.arctan2(*(end_row_col - ridge_row_col)))
        [curRow, curCol] = ridge_row_col

        ## start ridge follow
        # initialization    
        nCut = 48;
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
        ridge_cut = np.zeros([50, 41])
        ridge_pos = np.zeros([50, 2])
        ridge_pos[0, :] = np.array([curRow, curCol])
        ridge_angles = np.zeros(50)
        ridge_angles[0] = stepD
        ridge_dEdx = np.zeros(50)
        ridge_dE   = np.zeros(50)
        ridge_width = np.zeros(50)
        ridgeThreshold = 0.2/16. # threshold to stop ridge following
        cutThreshold = 0.2/16. # cut points below what value to assign 0
        for step_num in xrange(1, 50):
            cutAngle0 = adjust_angles(stepD+90) 
            cutAngles = np.arange(int(cutAngle0-nCut/2), cutAngle0+nCut/2+0.1)
            cutAngles = adjust_angles(cutAngles).astype('int')
            row_cuts = row_rotations[cutAngles] + curRow # row matrix 
            col_cuts = col_rotations[cutAngles] + curCol # col matrix
            cuts_shape = np.shape(col_cuts)
            energy_cuts = f_interp2d(row_cuts.flatten(), col_cuts.flatten(), grid=False).reshape(cuts_shape)
            energy_cuts[abs(energy_cuts)<cutThreshold] = 0
            dis = np.abs(col_0deg)
            width_cuts = np.sum(energy_cuts * dis, axis=1)/np.sum(energy_cuts, axis=1)
            w = np.min(width_cuts)  # min width
            index = np.argmin(width_cuts)   # min width index
            cutE = energy_cuts[index]      # cut energy vector
            cutEtot = np.sum(cutE)
            cutRow = row_cuts[index]
            cutCol = col_cuts[index]
            momentumRow = np.sum(cutRow*cutE)
            momentumCol = np.sum(cutCol*cutE)
            centerRow = momentumRow/cutEtot
            centerCol = momentumCol/cutEtot
            dEdx = cutEtot / (0.25 * 10.5) # in keV/um
            dE = f_interp2d(centerRow, centerCol, grid=False)
            cutD = cutAngles[index] # cut direction scalar
            stepD = (cutD-90)
            dRow = 0.25 * np.sin(deg2rad(stepD))
            dCol = 0.25 * np.cos(deg2rad(stepD))
            curRow = centerRow + dRow 
            curCol = centerCol + dCol

            ridge_pos[step_num] = np.array([curRow, curCol])
            ridge_angles[step_num] = stepD
            ridge_dEdx[step_num-1] = dEdx
            ridge_dE[step_num] = dE
            ridge_width[step_num-1] = w
            ridge_cut[step_num-1] = cutE

            if dE < ridgeThreshold:
                break
        ridge_pos = ridge_pos[:step_num+1][:-5] # ignore last N ridge points
        ridge_angles = ridge_angles[:step_num+1][:-5]
        ridge_dEdx = ridge_dEdx[:step_num]
        ridge_dE  = ridge_dE[:step_num+1]
        ridge_width = ridge_width[:step_num]
        ridge_cut = ridge_cut[:step_num]

        alpha_median = np.median(ridge_angles[1:]) 
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
        result['dE'] = ridge_dE
        result['width'] = ridge_width
        result['alpha_median'] = alpha_median
        result['alpha_linearReg'] = alpha_linearReg
        if int(file_idx) > 101:
            plotflag = False

        # plot
        if plotflag:

            fig = Figure()
            canvas = FigureCanvasAgg(fig)
            ax = fig.add_subplot(111)
            ax.imshow(image, interpolation='nearest', origin='lower')
            ax.scatter(ridge_pos[:,1], ridge_pos[:,0], c='g', edgecolor='none')
            ax.plot(true_pos['col'], true_pos['row'], 'ro-', markersize=1.5, markeredgewidth=0)
            ridge_col_min = np.floor(min(ridge_pos[:,1]))
            ridge_col_max = np.ceil(max(ridge_pos[:,1]))
            ridge_row_min = np.floor(min(ridge_pos[:,0]))
            ridge_row_max = np.ceil(max(ridge_pos[:,0]))
            ax.set_xlim([ridge_col_min-2.5, ridge_col_max+2.5])
            ax.set_ylim([ridge_row_min-2.5, ridge_row_max+2.5])
            ax.set_title('ridgeThreshold = %s || linear_reg: %.1f || median: %s  ' %(ridgeThreshold, alpha_linearReg, alpha_median))
            canvas.print_figure('%s/%s_%s.png' %(outdir, 'ridge', file_idx))

            fig = Figure()
            canvas = FigureCanvasAgg(fig)
            ax = fig.add_subplot(111)
            ax.plot(ridge_dE, 'bs-', lw=2)
            ax.set_xlabel('ridge #')
            ax.set_ylabel('keV')
            ax.set_title('dE', fontsize=24)
            canvas.print_figure('%s/%s_%s.png' %(outdir, 'dE', file_idx))

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
            ax.plot(np.arange(len(ridge_width)), ridge_width, 'bo-', lw=2)
            ax.set_xlabel('ridge #')
            ax.set_ylabel('pix')
            ax.set_title('Cut width')
            canvas.print_figure('%s/%s_%s.png' %(outdir, 'width', file_idx))
            
            fig = Figure()
            canvas = FigureCanvasAgg(fig)
            ax = fig.add_subplot(111)
            ax.plot(np.arange(np.shape(ridge_cut)[1]), np.transpose(ridge_cut))
            ax.set_xlabel('ridge #')
            ax.set_ylabel('keV')
            ax.set_title('cut energy')
            canvas.print_figure('%s/%s_%s.png' %(outdir, 'cutE', file_idx))
        
    
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
