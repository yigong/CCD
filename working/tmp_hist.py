figure()
hist(energy_list, bins=np.arange(4000, 200000, 100))
xlabel('ADC')
ylabel('counts')
ylim(0,300)
title('old median subtraction algorithm')

figure()
hist(energy_list_new, bins=np.arange(4000, 200000, 100))
xlabel('ADC')
ylabel('counts')
ylim(0,300)
title('new median subtraction algorithm')



show()


