import numpy as np, astropy, scipy, glob, sys
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit
from scipy.misc import factorial
import scipy.signal as sig
import scipy.stats as st
from scipy.interpolate import interp1d
import random
#import seaborn

import plot_conf_int
#import dosurvey, populate

#%pylab
#%matplotlib inline


def get_pop_pdfs(N_tots, alphas):
    # Function to get pdfs for population numbers

    pop_pdfs = np.full((len(alphas), len(N_tots)), 0.0)

    for ii, aa in enumerate(alphas):
        # Eq. 9 in the paper:
        pop_pdfs[ii] = aa ** 2 * N_tots * np.exp(-1 * aa * N_tots)

    return pop_pdfs


def get_rate_pdfs(N_tots, N_psr, alphas, t_life, f_b):
    # Function for calculating rates and rate pdfs

    # Initialize variables to hold the output
    rates = np.full((len(N_psr), len(N_tots)), 0.0)
    rate_pdfs = np.full(np.shape(rates), 0.0)

    for ii in range(len(N_psr)):
        # Eq. 10 in the paper:
        rates[ii] = N_tots * f_b[ii] / t_life[ii]

        # Eq. 11 in the paper:
        rate_pdfs[ii] = (alphas[ii] * t_life[ii] / f_b[ii]) ** 2 * rates[ii] * np.exp(
            -1 * (alphas[ii] * t_life[ii] / f_b[ii]) * rates[ii])

        # Normalize the above pdf:
        rate_pdfs[ii] = rate_pdfs[ii] / (scipy.integrate.simps(rate_pdfs[ii], rates[ii]))

    return rates, rate_pdfs


alphas_Pol = np.array([0.00668142, 0.01068674, 0.00652169, 0.00480014, 0.00701919, 0.00318349, 0.0091332, 0.01913401, 0.00743705])  #recr: Pol et al 2020
alphas_gen = np.array([0.00084191, 0.0033561, 0.00225884, 0.00186737, 0.00218487, 0.00122046, 0.00331682, 0.00350034, 0.00137995])  #with generic profile

N_psr_Pol = np.array([150, 94, 153, 208, 142, 314, 109, 52, 134])  #recr: Pol et al 2020
N_psr_gen = np.array([1188, 298, 443, 536, 458, 819, 301, 286, 725])  #generic profile

t_life = np.array([0.37e3, 2.93e3, 0.162e3, 0.293e3, 3.125e3, 0.24e3, 1.69e3, 0.06e3, 0.76e3])

f_b_Pol = np.array([5.72, 6.04, 4.59, 4.59, 4.59, 2.00, 4.59, 4.59, 4.59])
f_b_gen = np.array([1., 1., 1., 1., 1., 1., 1., 1., 1.])
f_b_guess = np.array([4, 4, 4, 4, 4, 4, 4, 6.9, 4])

N_tots = np.linspace(10, 4e3, 1000)




#Correct n_psr
#N_psr = np.array([281, 167, 214, 226, 205, 413, 165, 89])
#N_psr = np.array([53, 360])

#t_life = np.array([0.06e3, 0.06e3])


#f_b = np.array([fb_eff, 1.])

#pop_pdfs = get_pop_pdfs(N_tots, alphas)


#fb_eff=6.6

#f_b = np.array([fb_eff, 1.])

rates_Pol, rate_pdfs_Pol = get_rate_pdfs(N_tots, N_psr_Pol, alphas_Pol, t_life, f_b_Pol)
rates_gen, rate_pdfs_gen = get_rate_pdfs(N_tots, N_psr_gen, alphas_gen, t_life, f_b_gen)

#print(np.amax(rate_pdfs[0])-np.amax(rate_pdfs[1]))
'''
for ii, vals in enumerate(rate_pdfs_Pol):
    fb_eff = f_b_guess[ii]
    while np.abs(np.amax(rate_pdfs_Pol[ii])-np.amax(rate_pdfs_gen[ii]))>1e-3:
        fb_eff += 0.001
        f_b_gen[ii] = fb_eff
        rates_gen, rate_pdfs_gen = get_rate_pdfs(N_tots, N_psr_gen, alphas_gen, t_life, f_b_gen)
        print(fb_eff)
        names = np.array([r'h = 1, f$_b$ = {:.2f}'.format(fb_eff), r'spline, f$_b$ = 1'])
'''
ii=7

fb_eff = f_b_guess[ii]
while np.abs(np.amax(rate_pdfs_Pol[ii])-np.amax(rate_pdfs_gen[ii]))>1e-1:
    #print(np.abs(np.amax(rate_pdfs_Pol[ii])-np.amax(rate_pdfs_gen[ii])))
    fb_eff += 0.0001
    f_b_gen[ii] = fb_eff
    rates_gen, rate_pdfs_gen = get_rate_pdfs(N_tots, N_psr_gen, alphas_gen, t_life, f_b_gen)
    print(fb_eff)
    names = np.array([r'h = 1, f$_b$ = {:.2f}'.format(fb_eff), r'spline, f$_b$ = 1'])