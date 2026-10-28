# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 14:08:44 2026

@author: kgrun
"""


import os, glob, json
import numpy as np, scipy


import matplotlib.pyplot as pl
from scipy.optimize import curve_fit



def Poisson(x, lam):
    #return lam ** x * np.exp(-1 * lam) / scipy.special.factorial(x)
    return np.exp(np.log(lam)*x-lam-np.log(scipy.special.factorial(x)))

def line(x, a):
    return a * x


psr_names = ["B1534+12", "B1913+16",
             "J0509+3801", "J0737-3039A", "J1208-5936", "J1756-2251",
             "J1757-1854", "J1856-0039", "J1906+0746", "J1913+1102", 
             "J1946+2052"]


plotme = False

alpha_result = {}


runs_of_pop = 100
runs_per_pop = 100

n_bins = runs_of_pop + runs_per_pop      # total number of trials for the population
xaxis = np.arange(0, 200, 1)

result_dir = "./population_synthesis/popeval_res/"

for pname in psr_names:
    folderpath = result_dir + "pop_result_{}".format(pname)
    
    if os.path.isdir(folderpath) == False:
        print("No results for {} found.".format(pname))
        continue
    
    else:
        print("ANALYSE {}".format(pname))
    
    
    paths = sorted(glob.glob(folderpath + "/*npy"))

    # array holding the values for Ntot generated for each pulsar
    Ntot = np.arange(10, (len(paths) * 100 + 10), 100)



    # Preallocate arrays
    mean_prob = np.zeros((len(paths), n_bins))
    std_prob = np.zeros((len(paths), n_bins))
    mean_prob_norm = np.zeros((len(paths), n_bins))
    lambdas = np.zeros(len(paths))
    lambdas_err = np.zeros(len(paths))

    for jj, filename in enumerate(paths):
    
        # Shape is (runs_of_pop, runs_per_pop) = (100, 100)
        detect = np.load(filename).astype(int)
        n_detect = len(detect)  # runs_of_pop = 100
    
        prob = np.zeros((n_detect, n_bins))
    
        for ii, vals in enumerate(detect):
            counts = np.bincount(vals, minlength=n_bins)
            prob[ii] = counts / counts.sum()
    
        mean_prob[jj] = prob.mean(axis=0)
        std_prob[jj] = prob.std(axis=0)
    
        mean_prob_norm[jj] = mean_prob[jj] / mean_prob[jj].sum()
        
        popt, pconv = curve_fit(Poisson, xaxis, mean_prob_norm[jj], maxfev=10000)
    
        lambdas[jj] = popt[0]
        lambdas_fiterr = np.sqrt(pconv[0][0])

        lambdas_err[jj] = 1.96 * np.sqrt(popt[0] / Ntot[jj])

        '''
        if jj == 8:

            #pl.plot(xaxis, mean_prob[jj], color = 'r')
            pl.plot(xaxis, mean_prob_norm[jj], marker = 'o', color = 'k', label = 'Simulation')
            pl.plot(xaxis, Poisson(xaxis, *popt), marker = 's', color = 'c', label = 'Poisson')
            #pl.title(r"Results for $N_{\rm tot}$ = " + str(Ntot[jj]))
            pl.xlabel(r"N$_{\rm obs}$", fontsize = 22)
            pl.ylabel(r"P(N$_{\rm obs}$)", fontsize = 22)
            pl.legend(loc = 'best')
            pl.show()
            
        pl.errorbar(Ntot, lambdas, yerr = lambdas_err, marker = 'x', ls = '')

        pl.show()

        '''

    # fit the linear function to obtain alpha
    
    popt, pcov = curve_fit(line, Ntot, lambdas, sigma=lambdas_err)
    
    alpha_result[pname] = [popt[0], np.sqrt(pcov[0][0])]
    
    print(r"alpha = {} +/- {}".format(popt[0], np.sqrt(pcov[0][0])))
    print("")
    
    if plotme == True:
        pl.errorbar(Ntot, lambdas, yerr=lambdas_err, marker='x', ls='')
        pl.plot(Ntot, line(Ntot, *popt), color='g', label=r"$\alpha = $" + str(popt[0]))

        pl.xlabel(r"N$_{\rm tot}$")
        pl.ylabel(r"$\lambda$")

        pl.xlabel(r"Total number of DNS systems in Galaxy")
        pl.ylabel(r"Number of DNS systems detected")

        #pl.tight_layout()
        pl.legend(loc='best')
        pl.savefig(result_dir + 'alpha_{}.png'.format(pname))
        pl.show()


with open(result_dir + "alpha_dict.json", "w") as outfile:
    json.dump(alpha_result, outfile, indent=4)
    

outfile.close()
#N_pop = 1 / popt[0]

#print r"N_pop = {}".format(N_pop)
