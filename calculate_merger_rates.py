# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 11:52:33 2026

@author: kgrun
"""


import json
import numpy as np

#import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
 

import scipy
from scipy.interpolate import interp1d
from scipy.signal import fftconvolve


fsize=20

def calculate_ml_and_confint(x, y, conf_int = 0.95,
                             plot = False,
                             xlabel = r"$\cal{R}$$_{\rm MW}$ / Myr$^{-1}$", ylabel = r'P($\cal{R}$$_{\rm MW}$)',
                             color_fill = 'b', color_line = 'k', color_marker = 'b', label = '', display_conf_lines = True):
    
    #Script to plot the confidence intervals:
    #Calculate 95% limits from the posterior plot
    
    y_norm = y / scipy.integrate.simpson(y, x)
    total_area_under_post = scipy.integrate.simpson(y_norm, x)
    
    peak = x[np.argmax(y_norm)]
    
    #5 percent of this area:
    percent_5 = (1 - conf_int) * total_area_under_post / 2

    #Integrate starting from R = 0 till we hit 5 percent area and find that value of R
    for ii in range(1, len(y_norm) + 1):

        integral_ll = scipy.integrate.simpson(y_norm[: ii], x[: ii])

        if integral_ll >= percent_5:
            lower_limit_ind = ii

            break
            
    lower_limit = x[lower_limit_ind]
    
    reverse_post_R = y_norm[:: -1]
    reverse_possible_R = x[:: -1]

    for ii in range(1, len(reverse_post_R) + 1):

        integral_ul = scipy.integrate.simpson(reverse_post_R[: ii], x[: ii])

        if integral_ul >= percent_5:
            upper_limit_ind = len(y_norm) - ii

            break
    
    upper_limit = x[upper_limit_ind]
  
  
    if plot == True:
        plt.plot(x, y_norm, color = color_line, label = label)
    
        plt.fill_between(x[lower_limit_ind: upper_limit_ind], y_norm[lower_limit_ind: upper_limit_ind], color = color_fill, alpha = 0.2)
        if display_conf_lines:
            plt.axvline(x = lower_limit, ls = '--', color = color_marker, label = r"Lower limit = {:.2f}".format(lower_limit))
            plt.axvline(x = peak, ls = '--', color = color_marker, label = "Peak = {:.2f}".format(peak))
            plt.axvline(x = upper_limit, ls = '--', color = color_marker, label = "Upper limit = {:.2f}".format(upper_limit))

        plt.xlabel(xlabel, fontsize = 20)
        plt.ylabel(ylabel, fontsize = 20)

        plt.legend(loc = 'best')

        #plt.show()
    
    return peak, lower_limit, upper_limit




def get_pop_pdfs(N_tots, psrdict):
    """
    Function to get pdfs for population numbers
    """
    psr_keys = psrdict.keys()
    pop_pdfs = np.zeros( (len(psr_keys), len(N_tots)) )

    for ii, pkey in enumerate(psr_keys):
        # Eq. 9 in the paper:
        alpha = psrdict[pkey]["alpha"]
        pop_pdfs[ii] = alpha ** 2 * N_tots * np.exp(-1 * alpha * N_tots)

    return pop_pdfs


def get_rate_pdfs(N_tots, psrdict):
    """
    Function for calculating rates and rate pdfs
    """
    
    psr_keys = psrdict.keys()
    
    # Initialize variables to hold the output
    rates = np.zeros( (len(psr_keys), len(N_tots)) )
    rate_pdfs = np.zeros(np.shape(rates))

    for ii, pkey in enumerate(psr_keys):
        alpha  = psrdict[pkey]["alpha"]
        f_b    = psrdict[pkey]["f_b"]
        t_life = psrdict[pkey]["t_age"] + psrdict[pkey]["t_d"] 
        
        # Eq. 10 in the paper:
        rates[ii] = N_tots * f_b / t_life

        # Eq. 11 in the paper:
        rate_pdfs[ii] = (alpha * t_life / f_b) ** 2 * rates[ii] * np.exp( -1 * (alpha * t_life / f_b) * rates[ii])

        # Normalize the above pdf:
        rate_pdfs[ii] = rate_pdfs[ii] / (scipy.integrate.simpson(rate_pdfs[ii], rates[ii]))

    return rates, rate_pdfs


# recycled MSPs:    time for which the pulsar is detectable, t_d = merger time.
# J0737B & J1906:   non-recycled companion; t_d = time until PSR crosses death line
# 

# fb(1906) = 1 bc of profile implementation
# beaming correction factor for unknown fbs is # the average from the known fbs
pulsar_dict = {
    "B1534+12":    {"t_age": 208,   "t_d": 2700, "f_b": 6.0},
    "B1913+16":    {"t_age": 77,    "t_d": 301,  "f_b": 5.7},
    "J0509+3801":  {"t_age": 153,   "t_d": 576,  "f_b": 4.6},  # newly added by Pol et al. 2020
    "J0653+0443":  {"t_age": 6600,  "t_d": 2540, "f_b": 4.6},  # eccentric DNS from FAST
    "J0737-3039A": {"t_age": 159,   "t_d": 85,   "f_b": 2.0},
    "J1208-5936":  {"t_age": 2500,  "t_d": 7200, "f_b": 4.6},  # MMGPS discovery
    "J1756-2251":  {"t_age": 396,   "t_d": 1690, "f_b": 4.6},
    "J1757-1854":  {"t_age": 87,    "t_d": 76,   "f_b": 4.1},
    "J1846-0513":  {"t_age": 367,   "t_d": 8700, "f_b": 4.6},  # FAST
    "J1856-0039":  {"t_age": 586,   "t_d": 82,   "f_b": 4.6},  # FAST discovery
    "J1906+0746":  {"t_age": 0.11,  "t_d": 60,   "f_b": 1.0},  # modified pulse profile
    "J1913+1102":  {"t_age": 2625,  "t_d": 500,  "f_b": 4.6},
    "J1946+2052":  {"t_age": 247,   "t_d": 46,   "f_b": 2.1},
}



result_dir = "./population_synthesis/popeval_withGPPS/"
alphares_dict = json.load(open(result_dir + "alpha_dict.json", "r"))


for pulsar_name in pulsar_dict:
    if pulsar_name in alphares_dict:
        pulsar_dict[pulsar_name]["alpha"] = alphares_dict[pulsar_name][0]
    else:
        print(f"Warning: {pulsar_name} not found in alphares_dict")



## CALCULATE DISTRIBUTIONS #####

N_tots = np.linspace(1, 4e3, 2000)

pop_pdfs = get_pop_pdfs(N_tots, pulsar_dict)

rates, rate_pdfs = get_rate_pdfs(N_tots, pulsar_dict)







## PLOTTING #####



# Plot probability distributions for the individual populations

names = list(pulsar_dict.keys())

fig, ax = plt.subplots()
ax.set_prop_cycle(color = plt.rcParams['axes.prop_cycle'].by_key()['color'][:9],
                  linestyle = ['-', ':', '--', '-.', (0, (5, 5)), (0, (3, 1, 1, 1)), (0, (1, 5)), (0, (1, 1)), (0, (1, 3))])

for ii, vals in enumerate(pop_pdfs):
    ax.plot(N_tots, vals, label=names[ii])
    #plot_conf_int.plot_area_under_curve(N_tots, vals, xlabel = r"$N_{\rm det}$", ylabel = r"$P(N_{\rm det})$", label = names[ii])

ax.set_xscale('log')

ax.set_xlabel(r"$N_{\rm det}$", fontsize=fsize)
ax.set_ylabel(r"$P(N_{\rm det})$", fontsize=fsize)

ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
ax.grid(True, which='both')


plt.savefig(result_dir + "observable_pop_new.pdf", bbox_inches='tight')

plt.show()




'''
# Plot the TOTAL DNS distribution for each DNS system
for ii, vals in enumerate(pop_pdfs):
    #print names[ii]
    plt.plot(N_tots * f_b[ii], vals/f_b[ii], label=names[ii], ls=linestyles[ii])
    #plot_conf_int.plot_area_under_curve(N_tots * f_b[ii], vals/f_b[ii], xlabel = r"$N_{\rm tot}$", ylabel = r"$P(N_{\rm tot})$", label = names[ii])

plt.xscale('log')

plt.xlabel(r"$N_{\rm tot}$", fontsize=fsize)
plt.ylabel(r"$P(N_{\rm tot})$", fontsize=fsize)

plt.grid(True, which='both')
plt.tight_layout()

plt.legend(loc='best')

plt.savefig("total_population_new.pdf")

plt.show()
'''




###############################################################################
# PDFs for individual rate contributions                                      #
###############################################################################


fig, ax = plt.subplots()
ax.set_prop_cycle(color = plt.rcParams['axes.prop_cycle'].by_key()['color'][:9],
                  linestyle = ['-', ':', '--', '-.', (0, (5, 5)), (0, (3, 1, 1, 1)), (0, (1, 5)), (0, (1, 1)), (0, (1, 3))])


for ii, rpdf in enumerate(rate_pdfs):
    #print names[ii]
    ax.loglog(rates[ii], rpdf, label=names[ii])
    
    peak_rate, ll, ul = calculate_ml_and_confint(rates[ii], rpdf, 
                                                 plot = False,
                                                 label = names[ii], xlabel = r"$\cal{R}$", ylabel = r"P($\cal{R}$)")
    
    print("{:>12}: {:.2f} +{:.2f} -{:.2f}".format(names[ii], peak_rate, ul-peak_rate, peak_rate-ll))
    pulsar_dict[names[ii]]["rate"] = [peak_rate, ul-peak_rate, peak_rate-ll]



with open(result_dir + "pulsarval_dict.json", "w") as outfile:
    json.dump(pulsar_dict, outfile, indent=4)
outfile.close()
    


ax.set_ylim(1e-3, 50)
ax.set_xlim(1e-3, 50)

ax.set_xlabel(r"$\cal{R}$ / Myr$^{-1}$", fontsize=fsize)
ax.set_ylabel(r"P($\cal{R}$)", fontsize=fsize)

ax.grid(True, which='both')


ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")

plt.savefig(result_dir + "individual_rates_new.pdf", bbox_inches='tight')

plt.show()







# Take each DNS system; get integer rates; create rvs distribution:





def get_pmf_on_grid(rates_array, pdf_array, grid, interpolation="linear"):
    """Interpolate distribution onto grid and normalize."""

    if interpolation == "linear":
        x = rates_array
        x_grid = grid

    elif interpolation == "log":
        if np.any(rates_array <= 0) or np.any(grid <= 0):
            raise ValueError("Log interpolation requires positive rates.")

        x = np.log10(rates_array)
        x_grid = np.log10(grid)

    else:
        raise ValueError("interpolation must be 'linear' or 'log'.")

    func = interp1d(x, pdf_array, bounds_error=False, fill_value=0.0)
    
    pmf = np.clip(func(x_grid), 0, None)
    pmf_sum = pmf.sum()
    if pmf_sum <= 0:
        raise ValueError("Grid does not overlap with input distribution support.")

    return pmf / pmf_sum


def sum_all_rvs(rates_list, pdfs_list, dx=0.2, interpolation="linear"):

    lo = min(np.min(r) for r in rates_list)
    hi = max(np.max(r) for r in rates_list)

    n_points = int(round((hi - lo) / dx)) + 1
    common_grid = lo + dx * np.arange(n_points)

    pmfs = [ get_pmf_on_grid(r, p, common_grid, interpolation=interpolation) for r, p in zip(rates_list, pdfs_list)]
    pdf_final = pmfs[0]

    for pmf in pmfs[1:]:
        pdf_final = fftconvolve(pdf_final, pmf)
        pdf_final = np.clip(pdf_final, 0, None)
        pdf_final /= pdf_final.sum()

    rate_final = ( len(rates_list) * common_grid[0] + dx * np.arange(len(pdf_final))  )

    return rate_final, pdf_final


rates_final, pdf_final = sum_all_rvs(rates, rate_pdfs, dx=0.2)


calculate_ml_and_confint(rates_final, pdf_final, conf_int = 0.90, plot=True)
#plot_conf_int.plot_area_under_curve(rates_all, pdf_all)


   
plt.xscale('log')
plt.xlim(3, 100)

plt.grid(True, which = 'both')
plt.tight_layout()
plt.savefig(result_dir + "mw_merger_rate_new.pdf")
plt.show()




'''
### LIGO MERGER RATE EXTRAPOLATION ###
#Method from Koparappu, et al. 2008

d_r = 130    #Mpc
L_mw = 1.7
r_mw = rates_final / L_mw

ligo_rates = 7.4e-3 * (r_mw) * (2.26 * d_r / 100.) ** 3 #* 1.5

peak, lowlim, uplim = calculate_ml_and_confint(ligo_rates, pdf_final, conf_int = 0.90, xlabel = r'$\cal{R}$$_{\rm LIGO}$ / yr$^{-1}$', ylabel = r'P($\cal{R}$$_{\rm LIGO})$',
                                               plot=True)

plt.xscale('log')
plt.xlim(1, 40)

plt.grid(True, which = 'both')
plt.tight_layout()

plt.savefig(result_dir + "complete_LIGO_rate_new.pdf")
plt.show()




'''




'''
###############################################################################
#
# Comparison plot (similar to Fig. 8 in Grunthal et al.)
#
###############################################################################

#Make a nice 1-D plot to compare LIGO and our rate:
our_rates = np.array(peak)
ul = np.array([uplim])
ll = np.array([lowlim])

lf_rate = np.array([0.30]) #Luminosity function variation
lf_ll = np.array([0.10])
lf_ul = np.array([0.20])

ligo_rate = np.array([1.54]) #LIGO rate
ligo_ll = np.array([1.22])
ligo_ul = np.array([3.20])

bf_rate = np.array([0.35]) #Beaming factor variation, bf_new = 10
bf_ul = np.array([0.3])
bf_ll = np.array([0.13])

ellip_rate = np.array([0.31]) #Include elliptical galaxy contribution
ellip_ul = np.array([0.21])
ellip_ll = np.array([0.11])

allover_rate = np.array([0.45]) # All of the above corrections, except beaming factor, together
allover_ul = np.array([0.31])
allover_ll = np.array([0.16])

cli_rate = np.array([0.6]) #Chruslinska rate
cli_ul = np.array([1.2])
cli_ll = np.array([0.3])

map_rate = np.array([0.30]) #Mapelli rate
map_ul = np.array([0.30])
map_ll = np.array([0.28])


#Plot all the above
#plt.errorbar(np.arange(0, len(rates)), rates, yerr = [lol, upl], fmt = 'o', ls = '')
plt.figure(figsize = (20, 10))
l0a, caps1, c1 = plt.errorbar(our_rates, 0.01, xerr = [ll, ul], fmt = 'o', capsize = 6, label = 'This work')
l0b, caps2, c2 = plt.errorbar(ligo_rate, 0, xerr = [ligo_ll, ligo_ul], fmt = 'o', capsize = 6, label = 'LIGO DNS detection rate prediction')

l3, caps3, c3 = plt.errorbar(bf_rate, 0.1, xerr = [bf_ll, bf_ul], fmt = 'o', capsize = 6, label = 'Assume higher beaming correction factor')
l4, caps4, c4 = plt.errorbar(ellip_rate, 0.13, xerr = [ellip_ll, ellip_ul], fmt = 'o', capsize = 6, label = 'Including elliptical galaxy contribution')

l3, caps3, c3 = plt.errorbar(lf_rate, 0.15, xerr = [lf_ll, lf_ul], fmt = 'o', capsize = 6, label = 'Smaller mean luminosity')

l4, caps4, c4 = plt.errorbar(allover_rate, 0.2, xerr = [allover_ll, allover_ul], fmt = 'o', capsize = 6, label = 'Include elliptical correction and smaller mean luminosity')

plt.xscale('log')

plt.xlabel(r'$\cal{R}$ / yr$^{-1}$', fontsize = 40)

#This turns off the tick marks for y-axis
plt.setp(plt.gca(), yticks=())
plt.xticks(fontsize = 25)

plt.grid(True, which = 'both')

plt.legend(loc = 'best', frameon = True, fontsize = 16)
plt.tight_layout()

#plt.savefig("sim_var_new.pdf")

plt.show()


plt.figure(figsize = (20, 10))
l0a, caps1, c1 = plt.errorbar(our_rates, 0.01, xerr = [ll, ul], fmt = 'o', capsize = 6, label = 'This work')
l0b, caps2, c2 = plt.errorbar(ligo_rate, 0, xerr = [ligo_ll, ligo_ul], fmt = 'o', capsize = 6, label = 'LIGO DNS detection rate prediction')

plt.xscale('log')

plt.xlabel(r'$\cal{R}$ / yr$^{-1}$', fontsize = 40)

#Add the chrusylinski and mapelli merger rates:
l1, caps4, c4 = plt.errorbar(allover_rate, 0.1, xerr = [allover_ll, allover_ul], fmt = 'o', capsize = 6, label = 'Include elliptical correction and smaller mean luminosity')

l2 = plt.plot(0.0484, 0.13, marker = 's', markersize = 15, label = 'Reference model (Chruslinska et al., 2017)')
l3, caps4, c4 = plt.errorbar(cli_rate, 0.15, xerr = [cli_ll, cli_ul], fmt = 'o', capsize = 6, label = 'Most optimistic model (Chruslinska et al., 2017)')

l4, caps4, c4 = plt.errorbar(map_rate, 0.2, xerr = [map_ll, map_ul], fmt = 'o', capsize = 6, label = 'Mapelli & Giacobbo 2018')

#This turns off the tick marks for y-axis
plt.setp(plt.gca(), yticks=())
plt.xticks(fontsize = 25)

plt.grid(True, which = 'both')

plt.legend(loc = 'best', frameon = True, fontsize = 16)
plt.tight_layout()

#plt.savefig("comp_diff_rates_new.pdf")

plt.show()
'''
