edata = eData('./recon_tracks_45000.mat')
edata.betaMapping()
edata.alphaMapping()
edata.energyMapping()

BINS = linspace(-10, 10,21)

figure()
ax = gca()
counts, bins, poly_prev = hist(edata.back_projection, bins=BINS, histtype='step')
vlines(-3, 0, counts.max(), color='r',alpha=0.5, lw=1.5)
vlines(3, 0, counts.max(), color='r',alpha=0.5, lw=1.5)
xlabel('back projection at scattering window (mm)', fontsize='x-large')
ylabel('counts', fontsize='x-large')
ylim([0, 800])
title('No cut', fontsize='x-large')

lower_bound = edata.back_projection > -2.5
upper_bound = edata.back_projection < 2.5
in_window = np.logical_and(lower_bound, upper_bound).sum()
print '%s tracks are within window.'%(in_window)
print '%s tracks are left.'%(len(edata.diffusivity))


# end = 2
edata.cutData(edata.ends_num == 2)
figure()
ax = gca()
ax.plot(bins[1:], counts, drawstyle='steps', c='g', alpha=0.4, label='no cut')
vlines(-3, 0, counts.max(), color='r',alpha=0.5, lw=1.5)
vlines(3, 0, counts.max(), color='r',alpha=0.5, lw=1.5)
counts, bins, poly_prev = hist(edata.back_projection, bins=BINS, histtype='step',
                               color='b', label='num_ends=2')

xlabel('back projection position (mm)', fontsize='x-large')
ylabel('counts', fontsize='x-large')

title('number of ends cut', fontsize='x-large')
ylim([0, 800])
legend()
lower_bound = edata.back_projection > -2.5
upper_bound = edata.back_projection < 2.5
in_window = np.logical_and(lower_bound, upper_bound).sum()
print '%s tracks are within window.'%(in_window)
# diffusivity < 15
edata.cutData(edata.diffusivity < 15)
figure()
ax = gca()
ax.plot(bins[1:], counts, drawstyle='steps', c='g', alpha=0.4, label='num of ends cut')
vlines(-3, 0, counts.max(), color='r',alpha=0.5, lw=1.5)
vlines(3, 0, counts.max(), color='r',alpha=0.5, lw=1.5)
counts, bins, poly_tmp = ax.hist(edata.back_projection, bins=BINS, 
                                 histtype='step', color='b', label='diffusivity cut')
ax.legend()
xlabel('back projection position (mm)', fontsize='x-large')
ylabel('counts', fontsize='x-large')
title('diffusivity cut', fontsize='x-large')

alpha_flag = np.zeros_like(edata.diffusivity, dtype=bool)
for i in range(len(alpha_flag)):
    alpha = edata.alpha[i]
    energy = edata.track_energy[i]
    init_idx = edata.ends_idx[i]
    init_tmp = edata.ends_pos[i][init_idx]
    row_init = init_tmp[0]
    col_init = init_tmp[1]
    alpha_min = edata.alpha_min_map[row_init, col_init]
    alpha_max = edata.alpha_max_map[row_init, col_init]
    if alpha_min < alpha < alpha_max:
        alpha_flag[i] = True


# alpha cut
edata.cutData(alpha_flag)
figure()
ax = gca()
ax.plot(bins[1:], counts, drawstyle='steps', c='r', alpha=0.4, label='diffusivity cut')
vlines(-3, 0, counts.max(), color='r',alpha=0.5, lw=1.5)
vlines(3, 0, counts.max(), color='r',alpha=0.5, lw=1.5)
counts, bins, poly_tmp = ax.hist(edata.back_projection, bins=BINS, 
                                 histtype='step', color='b', label='direction cut')
ax.legend()
xlabel('back projection position (mm)', fontsize='x-large')
ylabel('counts', fontsize='x-large')
ylim([0, 280])
title('direction cut', fontsize='x-large')



# energy cut
energy_flag = np.zeros_like(edata.diffusivity, dtype=bool)
for i in range(len(energy_flag)):
    alpha = edata.alpha[i]
    energy = edata.track_energy[i]
    init_idx = edata.ends_idx[i]
    init_tmp = edata.ends_pos[i][init_idx]
    row_init = init_tmp[0]
    col_init = init_tmp[1]
    energy_min = edata.energy_min_map[row_init, col_init]
    energy_max = edata.energy_max_map[row_init, col_init]
    if energy_min < energy < energy_max:
        energy_flag[i] = True
        
edata.cutData(energy_flag)
figure()
ax = gca()
ax.plot(bins[1:], counts, drawstyle='steps', c='g', alpha=0.4, label='direction cut')
vlines(-3, 0, counts.max(), color='r',alpha=0.5, lw=1.5)
vlines(3, 0, counts.max(), color='r',alpha=0.5, lw=1.5)
counts, bins, poly_tmp = ax.hist(edata.back_projection, bins=BINS, 
                                 histtype='step', color='b', label='energy cut')

ax.legend()
xlabel('back projection position (mm)', fontsize='x-large')
ylabel('counts', fontsize='x-large')
ylim([0,50])
title('energy cut', fontsize='x-large')

lower_bound = edata.back_projection > -2.5
upper_bound = edata.back_projection < 2.5
in_window = np.logical_and(lower_bound, upper_bound).sum()

print '%s tracks are within window.'%(in_window)
show()

back_project_flag = np.zeros_like(edata.diffusivity, dtype=bool)
back_project_flag = np.logical_and(edata.back_projection>-2.5, edata.back_projection<2.5)
edata.cutData(back_project_flag)

close('all')

edata.pixelCoordinate()
edata.pointBack()
edata.measureAngles()










