from CCD.physics import Klein_Nishina_pdf
from CCD.physics import compton_gamma
from CCD.working.electron_hit import electron_hit



gamma_energy_keV = 1250
source_distance_cm = 7.5
activity_uCi = 150
time_min = 10
print '-'*40
print 'experiment setup:'
print '    gamma-ray energy: %s keV' %(gamma_energy_keV)
print '    source_distance: %s cm' %(source_distance_cm)
print '    source activity: %s uCi' %(activity_uCi)
print '    exposure time : %s min' %(time_min)

gamma_energy = gamma_energy_keV
distance = source_distance_cm * 10
# distance is in mm

activity = activity_uCi * 3.7e4
# activity is in bq. 3.7e4 is the conversion factor.

branching_ratio = 0.9998 * 2
time = 60 * time_min # sec

gamma_at_window = time*activity * branching_ratio* (2.5/distance)**2 *0.25
print '-'*40
print '- %.0f gammas reach at the mylar window.' % gamma_at_window

mylar_density = 1.38 # g/cm**3
mu = np.array([7.372e-2, 6.628e-2, 5.927e-2, 5.395e-2, 4.63e-2]) * mylar_density
energy = np.array([800, 1000, 1250, 1500, 2000])
mu_interp = np.interp(gamma_energy, energy, mu) # interaction per cm
mylar_thickness = 14e-4 # in cm
compton_rate = (1-np.exp(-mu_interp * mylar_thickness))

print '- the probability of compton scattering is %.5f.' % compton_rate

num_compton_tmp = gamma_at_window * compton_rate
num_compton = int(num_compton_tmp)

print '%.0f x %.5f = %s. \n- %s Compton electrons are produced.'  % (gamma_at_window, \
                                                            compton_rate, num_compton, num_compton)
        
if num_compton == 0:
    print '- the expected number of Compton electrons is %.3f.' % (num_compton_tmp)
    1/0

hit_probability = []

# num_compton = 100000
num_samples = 1000
pdf, sampled_cos_polar = Klein_Nishina_pdf(gamma_energy, num_samples)
cdf = pdf.cumsum()

geometry = {}
geometry['verticle_distance'] = 3
geometry['width'] = 36.8
geometry['length'] = 36.8
geometry['distance_from_window'] = 5

hit = 0
for i in range(num_compton):
    seed = np.random.random(1)
    tmp_bool = seed < cdf
#     print tmp_bool
#     print seed
    cos_polar = sampled_cos_polar[tmp_bool][0]
    compton_result = compton_gamma(gamma_energy, np.arccos(cos_polar))
    
    theta = compton_result['electron_scattered_angle']
    phi = np.random.random(1) * 2 * np.pi
    if electron_hit(theta, phi, geometry):
        hit += 1

average_hit_probability = 0.254
    
print 'we will get %s electron tracks on the CCD.' % (hit)
    

    
    
    
    