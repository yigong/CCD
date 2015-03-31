binsPB = np.arange(-20,20.5,0.5)
figPB, axPB = plt.subplots()

#axPB.hist(edataBeforeDiff.back_projection, bins=binsPB, histtype='stepfilled',\
#         align='left',color='g', label='Before width cut',alpha=0.5)
#axPB.hist(edataBeforeEnergy.back_projection, bins=binsPB, histtype='stepfilled',\
#         align='left',color='b', label='After width cut',alpha=0.5)
#axPB.hold(True)
axPB.hist(edataAfterEnergy.back_projection, bins=binsPB, histtype='stepfilled',\
         align='left',color='k',alpha=1)
axPB.legend()
axPB.set_xlim(-22, 22)
#axPB.set_ylim(0, 240)
axPB.set_xlabel('Pointing back position (mm)')
axPB.set_ylabel('Counts')

axPB.vlines([-2.75, 2.25], 0, 35, color='r',lw=1.5)
figPB.show()




ndim0 = ndim1 = 3506
rowMesh, colMesh = np.meshgrid(np.arange(ndim0), np.arange(ndim0), indexing='ij')
plMin = mlab.surf(rowMesh, colMesh, edata.energy_min_map)
plMax = mlab.surf(rowMesh, colMesh, edata.energy_max_map)
mlab.axes(xlabel='Row', ylabal='Column', zlabel='Energy')
mlab.outline(plMin)