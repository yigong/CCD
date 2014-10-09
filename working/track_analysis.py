import numpy as np
from scipy.io import loadmat
import pandas as pd
from CCD.plot import image2d


class DataCutter:
    ''' A class for cutting data arrays that are members of an inherited class
    
        An example usage is as follows:
        
        class Data(DataCutter):
            cutVars = ['var1', 'var2']
            
        Then calling:
        Data.cutVars( x>2 )
        where x is shape (N, ), the data will be cut according to the values of x
    '''
    
    # make a list of the data arrays to be cut by the cutData function 
    # these all need to be of shape (N, ...) where N is the number of data points
    
    cutVars = ['']
    
    def cutData(self, keep):
        ''' Cut the data using a logical keep array. This will cut all data
            arrays in the list self.cutVars
        '''
        
        print 'Cutting data from : %d to %d'%(keep.shape[0], keep.sum())
        for var in self.cutVars:
            if hasattr(self, var):
                ar = getattr(self, var)
                if np.shape(ar)[0] == keep.shape[0]:
                    ar = ar[keep]  
                    setattr(self, var, ar)     
                else:
                    print '%s not same shape: '%(var), ar.shape
            else:
                print '%s not found'%(var)
                
    def printData(self):
        ''' Print all the arrays which are cuttable '''
        
        print 'Variable name (shape):'
        for var in self.cutVars:
            if hasattr(self, var):
                ar = getattr(self, var)
                print '%s, '%(var, ), ar.shape
            else:
                print '%s not found'%(var)
                                
class eData(DataCutter):
    edge = 4

    def __init__(self, filename=None, numTrack=None):
        ''' Loading reconstructed electron tracks data '''
        if filename is None:
            print 'A new eData instance is created, but it\'s and empty instance.'
        else:
            
            recon_array = loadmat(filename)['recon_tracks'][0]
            dtype_names = recon_array.dtype.names
            recon_tracks = []
            if numTrack is None:
                numTrack = len(recon_array)
            print 'Loaded %s. There are %s tracks.'%(filename, numTrack)
            for i in range(numTrack): # loop the tracks
                recon_tmp = {} # the dict the save the track
                for j in range(len(dtype_names)): # loop the field
                    key = dtype_names[j]
                    recon_tmp[key] = np.squeeze(recon_array[i][key])
            
                lower_left = [np.min(recon_tmp['track_position'][0]), np.min(recon_tmp['track_position'][1])]
                recon_tmp['length'] = np.sum(recon_tmp['track_thinned'])
                recon_tmp['diffusivity'] = np.median(recon_tmp['ridge_FWHMum'])
                recon_tmp['ends_pos'] = recon_tmp['ends_pos'].astype(np.uint16) - 1 + lower_left# adjust for matlab 1 indexing
                recon_tmp['beta'] = recon_tmp['beta'].astype(np.float32)
                recon_tmp['ends_idx'] = recon_tmp['ends_idx'].astype(np.int16) - 1 # adjust for matlab 1 indexing
                recon_tmp['ends_num'] = recon_tmp['ends_num'].astype(np.int16)
            
                if recon_tmp['ends_num'] == 1:
                    recon_tmp['ends_pos'] = [recon_tmp['ends_pos']]
                    # to make sure recon_tmp['ends_pos']['ends_idx'] is a vector
                
                recon_tmp['px_num'] = np.sum(recon_tmp['track_1d'])
                
                end_idx = recon_tmp['ends_idx']
                # calculate back projected position on the window
                end_pos = recon_tmp['ends_pos'][end_idx]
                end_row = end_pos[0]
                end_col = end_pos[1]
                init_pos = [end_row, end_col]
                recon_tmp['init_pos'] = init_pos
                alpha = recon_tmp['alpha']
                
                if 0<=alpha<=90 or 270<=alpha<=360:
                    dim1 = (end_col - 1753 + 0.5) * 10.5/1000 + \
                            ((end_row+0.5)*10.5/1000 + self.edge) * \
                            np.tan(-1*alpha/180*np.pi)
                else:
                    dim1 = None    
                recon_tmp['back_projection'] = dim1
                recon_tracks.append(recon_tmp)
            dataFrame = pd.DataFrame(recon_tracks)
            self.cutVars = []
            for attr in dataFrame.columns:
                setattr(self, attr, dataFrame[attr].values)
                self.cutVars.append(attr)

    def betaMapping(self):
        from CCD.physics import compton_electron_array

        lower_left = [self.edge+0.5*10.5/1000, -1753.5*10.5/1000] # [dim0, dim1]
        px_size = 10.5/1000
        
        ndim0 = ndim1 = 3506
        dim0_vector, dim1_vector = np.meshgrid(np.arange(ndim0), np.arange(ndim0), indexing='ij')
        pos_window_coord = np.zeros([ndim0, ndim1])
        pos_window_coord_dim0 = lower_left[0] + dim0_vector*px_size
        pos_window_coord_dim1 = lower_left[1] + dim1_vector*px_size
        
        pos_window_coord = np.sqrt(pos_window_coord_dim0**2 +\
                                       pos_window_coord_dim1**2)
        
        
        self.beta_min_map = (np.arctan(1/pos_window_coord)) 
        self.beta_mean_map = (np.arctan(3.5/pos_window_coord))
        self.beta_max_map = (np.arctan(6/pos_window_coord))
        self.xVector = dim0_vector
        self.yVector = dim1_vector
        
    def alphaMapping(self):
        lower_left = [self.edge+0.5*10.5/1000, -1753.5*10.5/1000] # [dim0, dim1]
        px_size = 10.5/1000
        
        ndim0 = ndim1 = 3506
        dim0_vector, dim1_vector = np.meshgrid(np.arange(ndim0), np.arange(ndim0), indexing='ij')
        pos_window_coord = np.zeros([ndim0, ndim1])
        
        pos_window_coord_dim0 = lower_left[0] + dim0_vector*px_size
        pos_window_coord_dim1 = lower_left[1] + dim1_vector*px_size
        # alpha map calculated from geometry
        alpha_mean_map = np.arctan(-pos_window_coord_dim1/pos_window_coord_dim0)
        
        
        beta_sample = np.linspace(10, 40, 7)*np.pi/180.
        alpha_sigma_sample = np.array([9.80, 10.01, 10.52, 11.62, 12.96, 15.02, 18.15])
        
        alpha_sigma_map = np.interp(self.beta_mean_map.flatten(), beta_sample, alpha_sigma_sample).reshape((ndim0, ndim1))
        # for gaussian 95% --> 2 sigma each side
        self.alpha_min_map = alpha_mean_map - 2*alpha_sigma_map
        self.alpha_max_map = alpha_mean_map + 2*alpha_sigma_map
        self.alpha_mean_map = alpha_mean_map
        
    def energyMapping(self):
        from CCD.physics import compton_electron_array
        from CCD.misc.unit_convertion import radius_to_degree
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        
        lower_left = [self.edge+0.5*10.5/1000, -1753.5*10.5/1000] # [dim0, dim1]
        px_size = 10.5/1000
        
        ndim0 = ndim1 = 3506
        dim0_vector, dim1_vector = np.meshgrid(np.arange(ndim0), np.arange(ndim0), indexing='ij')
        pos_window_coord = np.zeros([ndim0, ndim1])
        pos_window_coord_dim0 = lower_left[0] + dim0_vector*px_size
        pos_window_coord_dim1 = lower_left[1] + dim1_vector*px_size
        
        # pos_window_coord = np.sqrt(pos_window_coord_dim0**2 +\
        #                                pos_window_coord_dim1**2)
        opposite = np.sqrt(3.5**2 + pos_window_coord_dim1**2)
        adjacent = pos_window_coord_dim0 
        self.phi_min_map = (np.arctan((opposite-2.5)/adjacent)) 
        self.phi_mean_map = (np.arctan(opposite/adjacent))
        self.phi_max_map = (np.arctan((opposite+2.5)/adjacent))
        
        self.energy_min_map = compton_electron_array(1173.2, self.phi_max_map.flatten()).reshape((ndim0, ndim1))
        self.energy_mean_map_1332 = compton_electron_array(1332.5, self.phi_mean_map.flatten()).reshape((ndim0, ndim1))
        self.energy_mean_map_1173 = compton_electron_array(1173.2, self.phi_mean_map.flatten()).reshape((ndim0, ndim1))
        self.energy_max_map = compton_electron_array(1332.5, self.phi_min_map.flatten()).reshape((ndim0, ndim1))
 
    def pixelCoordinate(self):
        ''' add attributes self.x/y/zEvent, which represent the position of the initial
            pixel in the detector coordinate system.
        '''
           
        pxSize = 10.5/1000
        x0 = self.edge + (0.5)*pxSize
        y0 = (-1753+0.5)*pxSize
        z0 = 0
        xEvent = []
        yEvent = []
        zEvent = []
        for [row, col] in self.init_pos:
            x = x0 + row*pxSize
            y = y0 + col*pxSize
            z = z0 + 0
            xEvent.append(x)
            yEvent.append(y)
            zEvent.append(z)
        self.xEvent = np.array(xEvent)
        self.yEvent = np.array(yEvent)
        self.zEvent = np.array(zEvent)
               
    def pointBack(self):
        ''' Use electron trajectory to point back where Compton scattering happening on the
            window. add attribute self.x/y/zWindow to the instance.
        '''
        
        self.pixelCoordinate()
        xWindow = []
        yWindow = []
        zWindow = []
        for (xEvent, yEvent, alpha) in zip(self.xEvent, self.yEvent, self.alpha):
            alphaRadius = alpha *np.pi/180.
            x = 0
            y = yEvent - xEvent*np.tan(alphaRadius)
            z = 3.5
            xWindow.append(x)
            yWindow.append(y)
            zWindow.append(z)
        self.xWindow = np.array(xWindow)
        self.yWindow = np.array(yWindow)
        self.zWindow = np.array(zWindow)
        
    def measureAngles(self):
        '''
        
        '''
        
        beta = []
        phi = []
        for (x, y, z) in zip(self.xWindow-self.xEvent, self.yWindow-self.yEvent, \
                             self.zWindow-self.zEvent):
            betaTmp = np.arctan(z/np.sqrt(x**2+y**2))
            phiTmp = np.arctan(np.sqrt(y**2+z**2)/np.abs(x))
            beta.append(betaTmp)
            phi.append(phiTmp)
        self.beta = beta
        self.phi = phi
        
    def energyInverting(self):
        Etrack = np.array(self.track_energy)
        Etrack = Etrack.astype(np.float32)
        m0C2 = np.ones_like(Etrack)*512
        cosPhi = np.cos(self.phi)
        sigmaC = np.sqrt(Etrack*(Etrack+2*m0C2))
        
        self.gammaEnergy = Etrack*m0C2/(cosPhi*sigmaC-Etrack)
    
    def on_mouse(self, event):
        event_row = int(event.ydata+0.5)
        event_col = int(event.xdata+0.5)
        
        iTrack = (self.image_segmented_filtered[event_row, event_col])
        info = "Pixel row,col:{}, {}\n".format(event_row, event_col)
        if iTrack != -1:
            #track_info = recon_tracks[segment_idx]
            alpha_radius = -1 * self.alpha[iTrack]/180.*np.pi
            end_idx = self.ends_idx[iTrack] 
            [end_row, end_col] = self.ends_pos[iTrack][end_idx]
             
            print self.image_filtered[event_row, event_col]
            # draw an arrow at the end 
            arrow_length = 10
            self.ax.annotate("", xy=(end_col+arrow_length*np.sin(alpha_radius), end_row-arrow_length*np.cos(alpha_radius)),
                        xytext=(end_col, end_row), xycoords='data', arrowprops=dict(arrowstyle='->',color='g' ) )
            # add info
            
            if self.back_projection[iTrack] != None:
                print_tmp = """backprojection: %.1f mm
alpha: %.1f degree
energy: %.1f keV
length: %s pixels
diffusivity: %.1f um
init segment: %s
num of ends : %s
index       : %s
""" % (self.back_projection[iTrack], self.alpha[iTrack],\
       self.track_energy[iTrack], self.length[iTrack],\
       self.diffusivity[iTrack], self.init_pos[iTrack],\
       self.ends_num[iTrack], iTrack)
            else:
                print_tmp = """backprojection: wrong direction
alpha: %.1f degree
energy: %.1f keV
length: %s pixels
diffusivity: %.1f um
init segment: %s
num of ends : %s
index       : %s
""" % (self.alpha[iTrack],\
       self.track_energy[iTrack], self.length[iTrack],\
       self.diffusivity[iTrack], self.init_pos[iTrack],\
       self.ends_num[iTrack], iTrack)
            
            info = info + print_tmp      
        else:
            #info = "Pixel row,col: %s, %s" % (event_row, event_col)
            print iTrack
            
        self.text.set_text(info)
        self.fig.canvas.draw()       

    def plotInteractive(self):
        
        image_filtered = np.zeros((3506, 3506))
        image_segmented_filtered = -1 * np.ones_like(image_filtered)
        image_segmented_filtered = image_segmented_filtered.astype(np.int16)
        for i, [rows, cols] in enumerate(self.track_position):
             image_filtered[rows, cols] += self.track_1d[i]
             image_segmented_filtered[rows, cols] = int(i)
        fig, ax = image2d(image_filtered)
        text = fig.text(0.8,0.5,'track_info')
        ax.set_title('filtered image')
        fig.show()
        self.image_filtered = image_filtered
        self.image_segmented_filtered = image_segmented_filtered
        self.text = text
        self.fig = fig
        self.ax = ax 
        
        fig.canvas.mpl_connect('button_press_event', self.on_mouse)
    
    def getTracks(self, index):
        ''' return a list of tracks. Each track is a eData instance '''
        
        from inspect import ismethod
        attrS = ['track_1d', 'track_energy', 'track_original', 'track_position','track_thinned']
        tracks = list()
        for i in index:
            oneTrack = eData()
            for attr in attrS:
                val = getattr(self, attr)
                setattr(oneTrack, attr, [val[i]])
            tracks.append(oneTrack)
        return tracks
    
    def 
        
        
        
        
        
        
        
        
        
        
        
        
