# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 17:16:09 2026

@author: kgrun
"""

import sys
sys.path.append("/u/kgrunthal/PsrPopPy2/lib/python") # ..
import numpy as np, astropy, scipy
import os
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit
from scipy.misc import factorial

import dosurvey as dosurvey, populate
import dosurvey_1906 as dosurvey_1906, populate_1906_spline


#Create a number of populations and see how many xx-like binaries we would detect


def combine_spline(beta,corr_func_MP,corr_func_IP):
    if beta < 10.:
        h = corr_func_MP(beta)
    
    else:
        h = max(corr_func_MP(beta),corr_func_IP(beta))
    
    return h




# Check if a pulsar name was provided
if len(sys.argv) < 2:
    print "Usage: python2 find_alpha.py <pulsar_name>"
    sys.exit(1)

# Get the pulsar name from command line argument
pulsar_name = sys.argv[1]



if pulsar_name == "J1906+0746":
    
    #new intensity profile correction
    profile_MP = np.loadtxt("./profile_files/mainpulse_mirror_log.txt").transpose()  #load determined profiles
    profile_IP = np.loadtxt("./profile_files/interpulse_mirror_log.txt").transpose()

    t_MP = np.append([-20.,-19,-18.,-17.,-16.,], np.append(profile_MP[0][36:58],[16.,17.,18.,19.,20.]))  #internal knots MP
    t_IP = np.append(profile_IP[0][2:6], np.append([-5.1,-5.,-4.,-3.,0.,3.,4.,5.,5.1], profile_IP[0][88:92])) + 18.82    #internal knots IP


    corr_func_MP = pol.LSQUnivariateSpline(profile_MP[0], profile_MP[1], t_MP)
    corr_func_IP = pol.LSQUnivariateSpline(profile_IP[0]+18.82, profile_IP[1], t_IP)
     
    corr_norm = max(corr_func_MP(np.arange(-12.,12.,0.001)))
      


pulsar_dict = {
    "B1534+12":    {"log_p": 1.5786392099680724,       "duty": 4.},
    "B1913+16":    {"log_p": 1.7709992051639407,       "duty": 16.9},
    "J0509+3801":  {"log_p": np.log10(76.5),           "duty": 18.},   # Pol et al. 2020
    "J0737-3039A": {"log_p": 1.3560258571931227,       "duty": 27.},
    "J0737-3039B": {"log_p": 3.4471580313422194,       "duty": None},
    "J1208-5936":  {"log_p": np.log10(28.7138598921),  "duty": 8.}, # Colom i Bernadich et al. 2023
    "J1756-2251":  {"log_p": np.log10(28.462),         "duty": 3.},
    "J1757-1854":  {"log_p": 1.3324384599156054,       "duty": 6.},
    "J1856-0039":  {"log_p": np.log10(23.39156058),    "duty": None}, # 2026
    "J1906+0746":  {"log_p": np.log10(144.073),        "duty": 1.},    # Grunthal et al. 2021
    "J1913+1102":  {"log_p": np.log10(27.28500685),    "duty": 6.},
    "J1946+2052":  {"log_p": 1.2278867046136734,       "duty": 6.},
}



surveys = ['PALFA_one_v_older', 'PHSURV', 'HTRU_low_1757', '1534_survey', 'PMSURV', 'GBNCC', "MMGPS-L", "GPPS"]

pops = np.arange(10, 2000, 100)
runs_of_pop = 100
runs_per_pop = 100

#lum_params = [-1.5, 0.94]     #From Bagchi, Lorimer, Jayanth, 2011
lum_params = [-1.1, 0.9]       #Fiducial parameters



model_name = './{}_z330.model'.format(pulsar_name)
#model_name = './1856_z330_varlum.model'

print "On {}".format(pulsar_name)

resfolder = "/scratch/kgrunthal/DNS-merger-rate/population_synthesis/test/"


detections = np.full((runs_of_pop, runs_per_pop), 0)

for xx, npop in enumerate(pops):

    filepath = resfolder + "/pop_result_{}/npop={:5d}.npy".format(pulsar_name, npop)  # check if the file already exists
    # print "On population:{}".format(npop)
    if os.path.exists(filepath):
        continue

    if pulsar_name == "J1906+0746":
        for ii in range(runs_of_pop):
        
            pop = populate_1906_spline.generate(npop, combine_spline, corr_func_MP, corr_func_IP, corr_norm, pDistPars = [log_p_1906_0746, 0.0], duty_percent = 1., orbits = False, nostdout = True, zscale = 0.33, siDistPars=[-1.4, 1.0], lumDistPars = lum_params)
            pop.write(outf = model_name)
        
        
            for jj in range(runs_per_pop):
                    
                population = dosurvey_1906.loadModel(model_name)
                    
                output = dosurvey_1906.run(population, pulsar_name, combine_spline, corr_func_MP, corr_func_IP, corr_norm, surveys, nostdout = True)
                    
                detections[ii][jj] = output[0][2].ndet + output[1][2].ndet + output[2][2].ndet + output[3][2].ndet + output[4][2].ndet + output[5][2].ndet
    
    
    else:
        for ii in range(runs_of_pop):
        
            pop = populate.generate(npop,
                                    pDistPars = [ pulsar_dict[pulsar_name]["log_p"], 0.0],
                                    duty_percent = pulsar_dict[pulsar_name]["duty"],
                                    orbits = False, nostdout = True,
                                    zscale = 0.33, siDistPars=[-1.4, 1.0],
                                    lumDistPars = lum_params)
            pop.write(outf = model_name)
        
            for jj in range(runs_per_pop):
                
                population = dosurvey.loadModel(model_name)
            
                output = dosurvey.run(population, pulsar_name, surveys, nostdout = True)
            
                detections[ii][jj] = output[0][2].ndet + output[1][2].ndet + output[2][2].ndet + output[3][2].ndet + output[4][2].ndet + output[5][2].ndet
            
    
    np.save(filepath, detections)
    
print "{} All Done!".format(pulsar_name)

