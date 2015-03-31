from copy import deepcopy
from matplotlib.pyplot import plotfile
import sys
plotFlag = False

edata = eData('../data/reconTracks_1st_single.mat')
edata.setSource([1173.2, 1332.5])
edata.betaMapping()
edata.alphaMapping()
edata.energyMapping()
#####
edataOriginal = deepcopy(edata)
#sys.exit("Error message")

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
edata.cutData(alpha_flag)
if plotFlag:
    titleInput = 'conservative cut on point back position'
    edata.positionHist(title=titleInput)
# energy inverting
edata.pixelCoordinate()
edata.pointBack()
#edata.measureAngles()
edata.energyInverting()
edata.addCutVar('gammaEnergy')
edata.addCutVar('phi')
#####
edataAfterPointBack = deepcopy(edata) # copy
#sys.exit("Error message")
if plotFlag:
    edata.electronEnergyPhiScatter(title=titleInput)
    edata.gammaEnergyPhiScatter(title=titleInput)

edata.cutData(edata.diffusivityMedian < 15)
#####
edataAfterDiff = deepcopy(edata)    # copy
#edataBeforeEnergy = deepcopy(edata) # copy
if plotFlag:
    titleInput = 'cut on diffusivity, number of ends, \
    point back position'
    edata.positionHist(title=titleInput)
    edata.electronEnergyPhiScatter(title=titleInput)
    edata.gammaEnergyPhiScatter(title=titleInput)


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
#####
edataAfterEnergy = deepcopy(edata) # copy
#edataBeforeEnds = deepcopy(edata)
if plotFlag:
    titleInput = 'cut on energy, diffusivity, number of ends, \
    point back position'
    edata.positionHist(title=titleInput)
    edata.electronEnergyPhiScatter(title=titleInput)
    edata.gammaEnergyPhiScatter(title=titleInput)

# edata.cutData(np.logical_or(edata.ends_num==2, edata.ends_num==3))
# if plotFlag:
#     titleInput = 'cut on number of ends, point back position'
#     edata.positionHist(title=titleInput)
#     edata.electronEnergyPhiScatter(title=titleInput)
#     edata.gammaEnergyPhiScatter(title=titleInput)
##edataAfterEnds = deepcopy(edata)

# back_project_flag = np.zeros_like(edata.diffusivity, dtype=bool)
# singleThickness = edata
# back_project_flag = np.logical_and(edata.back_projection>-1.5, edata.back_projection<1.5)
# edata.cutData(back_project_flag)
# titleInput = 'cut by scattering window width, energy, diffusivity,\
#  number of ends, point back position'
# edata.positionHist(title=titleInput)
# edata.electronEnergyPhiScatter(title=titleInput)
# edata.gammaEnergyPhiScatter(title=titleInput)
