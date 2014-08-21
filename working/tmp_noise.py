from numpy.fft import rfft


# sample rate gives us two things:
#    1. time interval: 1/sample rate
#    2. nyquest tells us that maximum in frequency axis is sample rate/2.
#  
# number of samples gives us:
#     1. length of time axis
#     2. step size in frequency axis
fft_flag = False

if fft_flag == True:
    sample_rate = 204868 # samples per second
    freq_max = sample_rate/2
    time_interval = 1./sample_rate
    
    num_samples = 3506
    t = np.arange(0, num_samples) * time_interval
    freq = np.linspace(0, freq_max, num_samples/2+1)
    
    
    result = np.zeros((3506, 1754))
    
    image = row_median_subtraction[1]
    fft_sig = rfft(image)
    fft_magnitude = np.abs(fft_sig)/1754

rows_to_use = np.arange(1, 1700, 300)

fig, (ax0, ax1) = plt.subplots(nrows=2)
ax0.plot(t, image[rows_to_use, :].T)
ax0.set_xlabel('Time(s)')
ax0.set_ylabel('ADC value')
ax0.set_title('Original Signal')
ax0.legend(['row # %s' % i for i in rows_to_use ])

ax1.plot(freq, fft_magnitude[rows_to_use, :].T)
ax1.set_xlabel('Frequency(Hz)')
ax1.set_ylabel('Intensity')
ax1.set_title('Frequency Spectrum')
ax1.legend(['row # %s' % i for i in rows_to_use ])
fig.suptitle('lower quadrant', fontsize=18)

rows_to_use = np.arange(1801, 3501, 300)

fig, (ax0, ax1) = plt.subplots(nrows=2)
ax0.plot(t, image[rows_to_use, :].T)
ax0.set_xlabel('Time(s)')
ax0.set_ylabel('ADC value')
ax0.set_title('Original Signal')
ax0.legend(['row # %s' % i for i in rows_to_use ])

ax1.plot(freq, fft_magnitude[rows_to_use, :].T)
ax1.set_xlabel('Frequency(Hz)')
ax1.set_ylabel('Intensity')
ax1.set_title('Frequency Spectrum')
ax1.legend(['row # %s' % i for i in rows_to_use ])
fig.suptitle('upper quadrant', fontsize=18)

plt.show()


# for i, row_data in enumerate(image):
#     if i%100 == 0:
#         print i
#     fft_output = rfft(row_data)/num_samples
#     magnitude_tmp =  [np.abs(val)/len(fft_output) for val in fft_output]
#     result[i] = magnitude_tmp
# 
# frequencies = [(k * 1.0/num_samples)*sample_rate for k in range(num_samples/2+1)]


    