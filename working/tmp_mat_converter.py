from CCD.io import save_mat
tmp = np.arange(12)
tmp= tmp.reshape(2,6)
tmp = [tmp]
save_mat('tmp', tmp, 'convert_test_2.mat')