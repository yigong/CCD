fig, ax = plt.subplots()
BINS = np.linspace(-10, 10, 41)
ax.hist(singleThickness.back_projection, bins=BINS, align='left', 
        label='14um', color='b', histtype='step', alpha=0.3)
ax.hist(doubleThickness.back_projection, bins=BINS+0.1, align='left',
        label='28um', color='r', histtype='step')
ax.legend()
ax.set_title('bin width = 0.5 mm')
plt.show()

# 
totalCounts_single = np.size(singleThickness.back_projection)
totalCounts_double = np.size(doubleThickness.back_projection)

single_3mm = np.sum(np.logical_and(singleThickness.back_projection < 1.5,
                                   singleThickness.back_projection > -1.5))

single_5mm = np.sum(np.logical_and(singleThickness.back_projection < 2.5,
                                   singleThickness.back_projection > -2.5))

double_3mm = np.sum(np.logical_and(doubleThickness.back_projection < 1.5,
                                   doubleThickness.back_projection > -1.5))

double_5mm = np.sum(np.logical_and(doubleThickness.back_projection < 2.5,
                                   doubleThickness.back_projection > -2.5))

print "14um: \n  # of events in the histogram(from -10 mm to 10 mm): %d"%(totalCounts_single)
print "  # of events within a diameter of 5 mm: %d"%(single_5mm)
print "  # of events within a diameter of 3 mm: %d"%(single_3mm)
print 50*'-'
print "28um: \n  # of events in the histogram(from -10 mm to 10 mm): %d"%(totalCounts_double)
print "  # of events within a diameter of 5 mm: %d"%(double_5mm)
print "  # of events within a diameter of 3 mm: %d"%(double_3mm)

