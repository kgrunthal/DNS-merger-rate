# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 11:14:35 2026

@author: kgrun
"""

import numpy as np
import numpy.random as nr

from scipy import stats

import matplotlib.pyplot as plt



def get_peak_and_hdi(samples, cred_interval=0.68):
    # Find the mode using KDE
    kde = stats.gaussian_kde(samples)
    x_grid = np.linspace(samples.min(), samples.max(), 2000)
    kde_vals = kde(x_grid)
    peak = x_grid[np.argmax(kde_vals)]
    
    # Sort samples for HDI calculation
    sorted_samples = np.sort(samples)
    n = len(sorted_samples)
    interval_idx = int(np.floor(cred_interval * n))
    
    # Find the shortest interval containing cred_interval fraction of samples
    n_intervals = n - interval_idx
    interval_widths = sorted_samples[interval_idx:] - sorted_samples[:n_intervals]
    min_idx = np.argmin(interval_widths)
    
    hdi_lower = sorted_samples[min_idx]
    hdi_upper = sorted_samples[min_idx + interval_idx]
    
    return peak, hdi_lower, hdi_upper



def Lambda(delta, i, Phi_SO):
    cosL = np.cos(delta)*np.cos(i) - np.sin(delta)*np.sin(i)*np.cos(Phi_SO)
    return np.arccos(cosL)



def Rho(alpha, beta, W10):
    t1 = np.cos(alpha)*np.cos(alpha+beta)
    t2 = np.sin(alpha)*np.sin(alpha+beta) * np.cos(W10/2.)
    return np.arccos(t1 + t2)


def fb_eff(alpha, rho):
    ul = np.array([min(alpha[i] +rho[i], np.pi/2) for i in range(len(alpha)) ])
    ll = np.array([max(0., alpha[i] - rho[i]) for i in range(len(alpha)) ])

    return np.abs((np.cos(ul) - np.cos(ll)))**(-1)





def fb_eff_1757(Ntrials):
    W10_draws = nr.uniform(18.0, 20.0, size= Ntrials) * (np.pi/180.)
    
    alpha1 = nr.normal(148.0, 3.0, size = Ntrials) * (np.pi/180.)
    alpha2 = nr.normal(144.0, 3.0, size = Ntrials) * (np.pi/180.)
    
    Phi_SO1 = nr.normal(168.0, 3.0, size = Ntrials) * (np.pi/180.)
    Phi_SO2 = nr.normal(168.0, 2.0, size = Ntrials) * (np.pi/180.)
    
    iota1 = 85.3 * (np.pi/180.)
    iota2 = 94.7 * (np.pi/180.)
    
    delta1 = nr.normal(46.0, 4.0, size = Ntrials) * (np.pi/180.)
    delta2 = nr.normal(52.0, 4.0, size = Ntrials) * (np.pi/180.)
    
    beta1 = np.pi - alpha1 - Lambda(delta1, iota1, Phi_SO1)
    beta2 = np.pi - alpha2 - Lambda(delta2, iota2, Phi_SO2)
    
    rho1 = Rho(alpha1, beta1, W10_draws)
    rho2 = Rho(alpha2, beta2, W10_draws)
    #print(rho1)
    
    fbeff1 = fb_eff(np.pi-alpha1, rho1)
    fbeff2 = fb_eff(np.pi-alpha2, rho2)

    
    return fbeff1, fbeff2




def fb_eff_1946(Ntrials):
    alpha = nr.normal(86.68, 0.11, size=Ntrials) * (np.pi/180.)
    rho_mp = nr.normal(23.39, 0.27, size=Ntrials) * (np.pi/180.)
    rho_ip = nr.normal(26.64, 0.4, size=Ntrials) * (np.pi/180.)
    
    fbeff_mp = fb_eff(alpha, rho_mp)
    fbeff_ip = fb_eff(alpha, rho_ip)
    
    
    return (fbeff_mp + fbeff_ip)/2



f1, f2 = fb_eff_1757(50000)


peak, lower, upper = get_peak_and_hdi(np.concatenate((f1,f2)))

print(f"Peak: {peak:.3f}")
print(f"-1σ: {peak - lower:.3f}  (lower bound: {lower:.3f})")
print(f"+1σ: {upper - peak:.3f}  (upper bound: {upper:.3f})")
print(f"Result: {peak:.3f} +{upper-peak:.3f}/-{peak-lower:.3f}")
print()
print()


fb1946 = fb_eff_1946(100000)

peak, lower, upper = get_peak_and_hdi(fb1946)

print(f"Peak: {peak:.3f}")
print(f"-1σ: {peak - lower:.3f}  (lower bound: {lower:.3f})")
print(f"+1σ: {upper - peak:.3f}  (upper bound: {upper:.3f})")
print(f"Result: {peak:.3f} +{upper-peak:.3f}/-{peak-lower:.3f}")




    