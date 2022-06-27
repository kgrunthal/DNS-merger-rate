# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 14:01:56 2020

@author: User
"""

import numpy as np, astropy, scipy, glob, sys
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as pl
#from scipy.optimize import curve_fit
#from scipy.misc import factorial
#import scipy.signal as sig
#import scipy.stats as st
#from scipy.interpolate import interp1d
#import random
#import seaborn

#import plot_conf_int



our_rates = np.array(4.19)
ul = np.array([2.88])
ll = np.array([1.32])

ligo_rate = np.array([10.11]) #LIGO rate
ligo_ll = np.array([4.48])
ligo_ul = np.array([15.59])

cli_rate = np.array([1.32]) #Chruslinska rate
cli_ul = np.array([1.32])
cli_ll = np.array([0.66])

#kru_rate = np.array([0.88]) #Kruckow rate #this was wrong, mail:16.08.2021
kru_rate = np.array([0.82])
kru_ul = np.array([3.68])
kru_ll = np.array([0.18])

ok_rate = np.array([3.20])
ok_ul = np.array([2.21])
ok_ll = np.array([0.88])

geom_rate = np.array([1.77])
geom_ul = np.array([2.76])
geom_ll = np.array([1.44])
#####################################
lf_rate = np.array([0.30]) #Luminosity function variation
lf_ll = np.array([0.10])
lf_ul = np.array([0.20])

ellip_rate = np.array([0.31]) #Include elliptical galaxy contribution
ellip_ul = np.array([0.21])
ellip_ll = np.array([0.11])

bf_rate = np.array([0.35]) #Beaming factor variation, bf_new = 10
bf_ul = np.array([0.3])
bf_ll = np.array([0.13])

allover_rate = np.array([0.45]) # All of the above corrections, except beaming factor, together
allover_ul = np.array([0.31])
allover_ll = np.array([0.16])

map_rate = np.array([0.30]) #Mapelli rate
map_ul = np.array([0.30])
map_ll = np.array([0.28])

#####################################

Pol_with_1906_rate = np.array([4.19])
Pol_with_1906_ul = np.array([2.98])
Pol_with_1906_ll = np.array([1.32])

Pol_recr_rate = np.array([3.86])
Pol_recr_ul = np.array([2.42])
Pol_recr_ll = np.array([1.21])

Pol_rate = np.array([4.17])
Pol_ul = np.array([2.63])
Pol_ll = np.array([1.3])


generic_rate = np.array([3.53])
generic_ul = np.array([2.10])
generic_ll = np.array([0.99])

ligo_rate = np.array([10.11]) #LIGO rate
ligo_ll = np.array([4.48])
ligo_ul = np.array([15.59])

#cli_rate = np.array([1.32]) #Chruslinska rate #this was wrong, mail:16.08.2021
#cli_ul = np.array([1.32])
#cli_ll = np.array([0.66])

cli_rate = np.array([5.52]) #Chruslinska rate #this was wrong, mail:16.08.2021
cli_ul = np.array([5.52])
cli_ll = np.array([2.76])

cli_base = np.array([0.45])

#kru_rate = np.array([0.88]) #Kruckow rate #this was wrong, mail:16.08.2021
kru_rate = np.array([0.82])
kru_ul = np.array([2.86]) #upper value: 3.68
kru_ll = np.array([0.64]) #lower value: 0.18




'''
#Plot all the above
#pl.errorbar(np.arange(0, len(rates)), rates, yerr = [lol, upl], fmt = 'o', ls = '')
#pl.figure(figsize = (20, 10))
l0a, caps1, c1 = pl.errorbar(our_rates, 0.01, xerr = [ll, ul], fmt = 'o', capsize = 6, label = 'This work')
l0b, caps2, c2 = pl.errorbar(ligo_rate, 0, xerr = [ligo_ll, ligo_ul], fmt = 'o', capsize = 6, label = 'LIGO DNS detection rate prediction')

l3, caps3, c3 = pl.errorbar(bf_rate, 0.1, xerr = [bf_ll, bf_ul], fmt = 'o', capsize = 6, label = 'Assume higher beaming correction factor')
l4, caps4, c4 = pl.errorbar(ellip_rate, 0.13, xerr = [ellip_ll, ellip_ul], fmt = 'o', capsize = 6, label = 'Including elliptical galaxy contribution')

l3, caps3, c3 = pl.errorbar(lf_rate, 0.15, xerr = [lf_ll, lf_ul], fmt = 'o', capsize = 6, label = 'Smaller mean luminosity')

l4, caps4, c4 = pl.errorbar(allover_rate, 0.2, xerr = [allover_ll, allover_ul], fmt = 'o', capsize = 6, label = 'Include elliptical correction and smaller mean luminosity')

pl.xscale('log')

pl.xlabel(r'$\cal{R}$ [yr$^{-1}$]', fontsize = 40)

#This turns off the tick marks for y-axis
pl.setp(pl.gca(), yticks=())
pl.xticks(fontsize = 25)

pl.grid(True, which = 'both')

pl.legend(loc = 'best', frameon = True, fontsize = 16)
pl.tight_layout()

pl.savefig("sim_var_geom.pdf")

pl.show()
'''

pl.figure(figsize = (15, 8))
#pl.figure(figsize = (500/300, 500/300))

l3, caps4, c4 = pl.errorbar(cli_rate, 0.04, xerr = [cli_ll, cli_ul], fmt = 'o', capsize = 6, label = 'most optimistic model by Chruslinska et al., 2017')
l6, caps6, c6 = pl.errorbar(cli_base, 0.04, marker='s', ls='', markersize=8, label = 'reference model by Chruslinska et al., 2017')
l5, caps5, c5 = pl.errorbar(kru_rate, 0.035, xerr = [kru_ll, kru_ul], fmt = 'o', capsize = 6, label = 'models by Kruckow et al., 2018')

l1a, caps4, c4 = pl.errorbar(Pol_recr_rate, 0.02, xerr = [Pol_recr_ll, Pol_recr_ul], fmt = 'o', capsize = 6, label = 'Pol et al. 2020 scheme')
l0a, caps1, c1 = pl.errorbar(Pol_with_1906_rate, 0.015, xerr = [Pol_with_1906_ll, Pol_with_1906_ul], fmt = 'o', capsize = 6, label = 'enhancing J1906+0746')
l4, caps4, c4 = pl.errorbar(generic_rate, 0.01, xerr = [generic_ll, generic_ul], fmt = 'o', capsize = 6, label = 'generic profile')

l0b, caps2, c2 = pl.errorbar(ligo_rate, 0, xerr = [ligo_ll, ligo_ul], fmt = 'o', capsize = 6, label = 'LIGO DNS merger detection rate prediction')

#l1b, caps4, c4 = pl.errorbar(Pol_rate, 0.03, xerr = [Pol_ll, Pol_ul], fmt = 'o', capsize = 6, label = 'Pol et al. 2020')

#This turns off the tick marks for y-axis
pl.setp(pl.gca(), yticks=())

pl.xlabel(r'$\cal{R}$$_\mathrm{LIGO}$ / yr$^{-1}$', fontsize = 40)
pl.xticks(fontsize = 25)
pl.xscale('log')


pl.grid(True, which = 'both')

pl.legend(loc = 'best', frameon = True, fontsize = 16)
pl.tight_layout()

pl.savefig("DNS_merger_rate_compilation.pdf")

pl.savefig("DNS_merger_rate_compilation.png", bbox_inches='tight', dpi=300)
pl.show()
