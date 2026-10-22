import sys
sys.path.append("../../PsrPopPy2-master/lib/python")
import numpy as np, astropy, scipy
import os
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit
from scipy.misc import factorial
import scipy.interpolate as pol

import dosurvey_1906_stepfunc as dosurvey, populate_1906

#Create a number of populations and see how many 1757-like binaries we would detect
def main():
    
    rho = float(sys.argv[1])
    
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
    #profile = np.loadtxt("mainpulse_mirror_log.txt").transpose()        
    #corr_func = pol.interp1d(profile[0][15:79],profile[1][15:79],kind='cubic', bounds_error=False, fill_value= -1.2)
    #corr_norm = max(corr_func(np.arange(-12.,12.,0.1)))
    
    
    #lum_params = [-1.5, 0.94]     #From Bagchi, Lorimer, Jayanth, 2011
    lum_params = [-1.1, 0.9]       #Fiducial parameters
    
    detections = np.full((runs_of_pop, runs_per_pop), 0)
    
    model_name = './1906_0746_z330_{}.model'.format(int(rho))
    #model_name = './1906_0746_z330_varlum.model'
    print "On 1906+0746"
    
        
    for xx, npop in enumerate(pops):
    
    
        filepath = "./1906_step_func_{}/npop={:5d}.npy".format(int(rho), npop)  # check if the file already exists
    
        # print "On population:{}".format(npop)
        if os.path.exists(filepath):
            print "Pfad {}.npy existiert".format(npop)
            continue
    
    
        for ii in range(runs_of_pop):
    
            pop = populate_1906.generate(npop, rho, pDistPars = [log_p_1906_0746, 0.0], duty_percent = 1., orbits = False, nostdout = True, zscale = 0.33, siDistPars=[-1.4, 1.0], lumDistPars = lum_params)
            pop.write(outf = model_name)
    
    
            for jj in range(runs_per_pop):
                
                population = dosurvey.loadModel(model_name)
                
                output = dosurvey.run(population, rho, surveys, nostdout = True)
                
                detections[ii][jj] = output[0][2].ndet + output[1][2].ndet + output[2][2].ndet + output[3][2].ndet + output[4][2].ndet + output[5][2].ndet
                
        #filepath = "./pop_result_z330_1906_0746/npop={:5d}.npy".format(npop)
        np.save(filepath, detections)
        
    print "1906+0746 All Done!"




if __name__ == '__main__':
    main()
