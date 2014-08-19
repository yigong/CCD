from numpy.fft import rfft

sample_rate = 204868 # samples per second
num_samples = 3506

result = np.zeros((3506, 1754))

image = row_median_subtraction[1]

for i, row_data in enumerate(image):
    print i 
    fft_output = rfft(row_data)
    magnitude_tmp =  [np.abs(val)/len(fft_output) for val in fft_output]
    result[i] = magnitude_tmp

frequencies = [(k * 1.0/num_samples)*sample_rate for k in range(num_samples/2+1)]
    