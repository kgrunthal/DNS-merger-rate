import numpy as np, astropy, scipy, glob, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit
from scipy.misc import factorial
import scipy
#import scipy.signal as sig
#import scipy.stats as st
#from scipy.interpolate import interp1d
#import seaborn

#import plot_conf_int
#import dosurvey, populate

#%pylab
#%matplotlib inline


def Poisson(x, lam):
    #return lam ** x * np.exp(-1 * lam) / scipy.special.factorial(x)
    return np.exp(np.log(lam)*x-lam-np.log(scipy.special.factorial(x)))


paths = sorted(glob.glob("./pop_result_1906_0746_random_alpha/*npy"))

Ntot = np.arange(10, (len(paths) * 100 + 10), 100)
xaxis = np.arange(0, 200, 1)

detect = np.load(paths[0])
prob = np.full((len(detect), 200), 0.0)
mean_prob = np.full((len(paths), 200), 0.0)
std_prob = np.full((len(paths), 200), 0.0)
mean_prob_norm = np.full((len(paths), 200), 0.0)

lambdas = np.full((len(paths)), 0.0)
lambdas_err = np.full((len(paths)), 0.0)

for jj, filename in enumerate(paths):

    prob = np.full((len(detect), 200), 0.0)

    print "On filename:{}".format(filename)
    detect = np.load(filename)
    detect = np.asarray(detect, dtype=np.int)

    for ii, vals in enumerate(detect):
        unique, counts = np.unique(vals, return_counts=True)

        probability = counts / np.sum(counts, dtype=np.float64)

        prob[ii][unique] = probability

    mean_prob[jj] = np.mean(prob, axis=0)
    std_prob[jj] = np.std(prob, axis=0)

    mean_prob_norm[jj] = mean_prob[jj] / np.sum(mean_prob[jj])
    #print mean_prob_norm[jj]

    popt, pconv = curve_fit(Poisson, xaxis, mean_prob_norm[jj], maxfev=10000)
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
      '''
    lambdas[jj] = popt[0]
    lambdas_err[jj] = 1.96 * np.sqrt(popt[0] / len(Ntot))

pl.errorbar(Ntot, lambdas, yerr = lambdas_err, marker = 'x', ls = '')

pl.show()

#print lambdas

#print r"err = {}".format(lambdas_err)

def line(x, a):
    return a * x


popt, pcov = curve_fit(line, Ntot, lambdas, sigma=lambdas_err)

print r"$\alpha$ = {}".format(popt)
print r"alphaerr = {}".format(np.sqrt(pcov[0][0]))

pl.errorbar(Ntot, lambdas, yerr=lambdas_err, marker='x', ls='')
pl.plot(Ntot, line(Ntot, *popt), color='g', label=r"$\alpha = $" + str(popt[0]))

pl.xlabel(r"N$_{\rm tot}$")
pl.ylabel(r"$\lambda$")

pl.xlabel(r"Total number of DNS systems in Galaxy")
pl.ylabel(r"Number of DNS systems detected")

#pl.tight_layout()
pl.legend(loc='best')
pl.savefig('alpha_1906_0746_random_alpha.png')
pl.show()

N_pop = 1 / popt[0]

print r"N_pop = {}".format(N_pop)
