from CCD.segmentation import noise_sigma
from CCD.io import load_fits
from CCD.segmentation import segment_tracks
from CCD.io import save_mat
base_dir = '/Users/Yigong/Research/photon_diagnostics/data/08_18_14_Am241/measurement'
tmp_list = []
for file_idx in [3]:

    GAIN = [5.6090e-4, 5.4834e-4]
    image_ADC = load_fits(base_dir + '/image_%s_corrected.fit' % file_idx)
    print 'Segmenting file: %simage_%s_corrected.fit\n' % (base_dir, file_idx)
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
    save_mat('segmented_tracks_%s' % file_idx, tmp, '%ssegmented_tracks_%s' % (base_dir, file_idx))

save_mat('segmented_tracks_%s' % file_idx, tmp_list, '%ssegmented_tracks_%s' % (base_dir, file_idx))    