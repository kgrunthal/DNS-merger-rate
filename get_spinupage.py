# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 09:50:26 2026

@author: kgrun
"""

import numpy as np




def birthperiod(pdict, alpha = 1.1e-15):
    p0 = (pdict["P"]*pdict["Pdot"] / alpha)**(3/7)
    return p0


def spinupline_age(pdict, pbirth):
    fac = (pbirth/pdict["P"])**2
    t_age = pdict["tau_c"] * (1-fac)
    return t_age
    


psr_dict = {"J0653+0443": {"P": 0.0475363055452, "Pdot": 1.15e-19, "tau_c": 6600},
            "J1846-0513": {"P": 0.0233695392907, "Pdot": 1.0106e-18, "tau_c":(10**8.56)*1e-6},
            "J1856-0039": {"P": 0.0233915605837, "Pdot": 6.532e-19, "tau_c": 568},
    }




for k in psr_dict.keys():
    p0 = birthperiod(psr_dict[k])
    t = spinupline_age(psr_dict[k], p0)
    
    print(k, np.round(t), np.round(psr_dict[k]["tau_c"]))