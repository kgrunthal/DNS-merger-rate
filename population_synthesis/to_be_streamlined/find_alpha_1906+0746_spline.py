import sys
sys.path.append("../../PsrPopPy2-master/lib/python")
import numpy as np, astropy, scipy
import os
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit
from scipy.misc import factorial
import scipy.interpolate as pol

import dosurvey_1906 as dosurvey, populate_1906_spline

#Create a number of populations and see how many 1757-like binaries we would detect

    
log_p_1102 = np.log10(27.28500685)
log_p_1946 = 1.2278867046136734
log_p_0737A = 1.3560258571931227
log_p_1757 = 1.3324384599156054
log_p_1534 = 1.5786392099680724
log_p_1913 = 1.7709992051639407
log_p_0737B = 3.4471580313422194
log_p_1756_2251 = np.log10(28.462)
log_p_1906_0746 = np.log10(144.073)
    
#survey_dict = {'1946': 'PALFA_one_v_older', '0737A': 'PHSURV', '1757': 'HTRU_low_1757', '1534': '1534_survey', '1913': '1913_survey', '0737B': 'PHSURV', '1102': 'PALFA_one_v_older'}

surveys = ['PALFA_one_v_older', 'PHSURV', 'HTRU_low_1757', '1534_survey', 'PMSURV', 'GBNCC']
pops = np.arange(10, 2000, 100)   #default: 10,6000,100
runs_of_pop = 100  # type: int, default 100
runs_per_pop = 100
    


#new intensity profile correction
profile_MP = np.loadtxt("mainpulse_mirror_log.txt").transpose()  #load determined profiles
profile_IP = np.loadtxt("interpulse_mirror_log.txt").transpose()

t_MP = np.append([-20.,-19,-18.,-17.,-16.,], np.append(profile_MP[0][36:58],[16.,17.,18.,19.,20.]))  #internal knots MP
t_IP = np.append(profile_IP[0][2:6], np.append([-5.1,-5.,-4.,-3.,0.,3.,4.,5.,5.1], profile_IP[0][88:92])) + 18.82    #internal knots IP


corr_func_MP = pol.LSQUnivariateSpline(profile_MP[0], profile_MP[1], t_MP)
corr_func_IP = pol.LSQUnivariateSpline(profile_IP[0]+18.82, profile_IP[1], t_IP)
 
corr_norm = max(corr_func_MP(np.arange(-12.,12.,0.001)))
  
def combine_spline(beta,corr_func_MP,corr_func_IP):
    if beta < 10.:
        h = corr_func_MP(beta)
    
    else:
        h = max(corr_func_MP(beta),corr_func_IP(beta))
    
    return h



#lum_params = [-1.5, 0.94]     #From Bagchi, Lorimer, Jayanth, 2011
lum_params = [-1.1, 0.9]       #Fiducial parameters
    
detections = np.full((runs_of_pop, runs_per_pop), 0)

model_name = './1906_0746_spline.model'
#model_name = './1906_0746_z330_varlum.model'
print "On 1906+0746"
    
        
for xx, npop in enumerate(pops):
    
    filepath = "./pop_result_1906_0746_new/npop={:5d}.npy".format(npop)  # check if the file already exists
    
    # print "On population:{}".format(npop)
    if os.path.exists(filepath):
        print "Pfad {}.npy existiert".format(npop)
        continue
    
    for ii in range(runs_of_pop):
    
        pop = populate_1906_spline.generate(npop, combine_spline, corr_func_MP, corr_func_IP, corr_norm, pDistPars = [log_p_1906_0746, 0.0], duty_percent = 1., orbits = False, nostdout = True, zscale = 0.33, siDistPars=[-1.4, 1.0], lumDistPars = lum_params)
        pop.write(outf = model_name)
    
    
        for jj in range(runs_per_pop):
                
            population = dosurvey.loadModel(model_name)
                
            output = dosurvey.run(population, combine_spline, corr_func_MP, corr_func_IP, corr_norm, surveys, nostdout = True)
                
            detections[ii][jj] = output[0][2].ndet + output[1][2].ndet + output[2][2].ndet + output[3][2].ndet + output[4][2].ndet + output[5][2].ndet
                
    #filepath = "./pop_result_z330_1906_0746/npop={:5d}.npy".format(npop)
    np.save(filepath, detections)
    
print "1906+0746 All Done!"

