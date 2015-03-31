from CCD.segmentation import noise_sigma
from CCD.io import load_fits
from CCD.segmentation import segment_tracks
from CCD.io import save_mat
from CCD.calibrate import self_subtraction

base_dir = '/Users/Yigong/Research/photon_diagnostics/data/10_28_14_1secFarAway/measurement'
tmp_list = []
GAIN = [5.6090e-4, 5.4834e-4]

for file_idx in range(1,74):
    file_name = base_dir + '/image_%s.fit' % file_idx
    image_ADC = self_subtraction(file_name, save_flag=False)
    
    print 'Segmenting file: %s/image_%s.fit\n' % (base_dir, file_idx)
    nrow, ncol = np.shape(image_ADC)
    
    image_bottom_ADC = image_ADC[:nrow/2]
    image_top_ADC = image_ADC[nrow/2:]
    
    sigma_bottom_ADC = noise_sigma(image_bottom_ADC)['STD']
    sigma_top_ADC = noise_sigma(image_top_ADC)['STD']
    
    image = np.concatenate((image_bottom_ADC * GAIN[0], image_top_ADC * GAIN[1]))
    
    sigma_bottom = sigma_bottom_ADC * GAIN[0]
    sigma_top = sigma_top_ADC * GAIN[1]
    THRESHOLDS = [10*sigma_bottom, 10*sigma_top]
    
    tmp = segment_tracks(image, THRESHOLDS)
    tmp_list.extend(tmp)
    save_mat('segmented_tracks_%s' % file_idx, tmp, '%s/segmented_tracks_%s' % (base_dir, file_idx))

save_mat('segmented_tracks_all', tmp_list, '%s/segmented_tracks_farAway'%(base_dir))    