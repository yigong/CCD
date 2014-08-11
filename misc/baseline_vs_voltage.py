import numpy as np
import matplotlib.pyplot as plt
voltages = np.arange(10,150,10)
L1 = [65535, 65535, 6588, 3975, 4801, 5509, 6293, 6971, 7617, 8309, 8958, 9680, 10335, 10916]
L2 = [65535, 65535, 5597, 5018, 5742, 6213, 6703, 7101, 7500, 7917, 8312, 8736, 9143, 9523]
U1 = [26592, 27616, 27616, 27616, 27616, 27616, 27616, 27616, 27616, 27520, 27360, 27136, 27136, 27136]
U2 = [65535, 5777, 4977, 4096, 4096, 2442, 2251, 2199, 2092, 2116, 2097, 2307, 2449, 2629]
data_160 = [L1, L2, U1, U2]

L1 = [65535, 65535, 8489, 8516, 4973, 5708, 6420, 7085, 7737, 8435, 9097, 9813, 10463, 11040]
L2 = [65535, 65535, 5822, 5864, 5702, 6240, 6692, 7074, 7484, 7913, 8321, 8733, 9169, 9549]
U1 = [26592, 27616, 27616, 27616, 27616, 27136, 27392, 27584, 27392, 27136, 27104, 27104, 27104, 27104]
U2 = [65535, 7281, 6321, 6289, 4433, 3095, 2730, 2560, 2286, 2208, 2217, 2291, 2436, 2616]
data_200 = [L1, L2, U1, U2]

quad_str = ['L1', 'L2', 'U1', 'U2']
fig, ax = plt.subplots(2,2)
for i, a in enumerate(ax.flatten('F')):
    a.plot(voltages, data_160[i], 'bx-',voltages, data_200[i], 'rx-')
    a.legend(["160K", "200K"])
    a.set_title(quad_str[i])
    a.set_xlabel('Bias Voltage[V]')
    a.set_ylabel('ADC value')
    a.grid()
    a.set_ylim((0, 66000))

plt.show()