
from CCD.io.fits import save_fits, load_fits


image_full = load_fits('/Users/Yigong/Documents/beamCCD/data/05_02_14_images_to_show/beta_beam_window.fit')
beta_image = image_full[1766:3530, 1764:3528]

left = beta_image[..., :882]
right = beta_image[..., 882:]

flipped = np.concatenate((right, left), 1).astype(np.int32)

#save_fits(beta_image, '')
#save_fits(flipped, '')
