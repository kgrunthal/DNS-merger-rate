[![DOI](https://zenodo.org/badge/147532251.svg)](https://zenodo.org/badge/latestdoi/147532251)

# On the detectability of eccentric binary pulsars, Bagchi, Lorimer, Wolfe, 2013, MNRAS

This repository holds the Fortran code written and developed by Bagchi, et al. in their analysis on the degradation factor for pulsars in binary orbits ([link to paper](http://adsabs.harvard.edu/abs/2013MNRAS.432.1303B)). The original host for the code on the WVU webpage is not available, and hence the code has been moved here with permission from the authors. This code is publicly available for use contingent on appropriate credit and citation to the paper above.

The instructions on using the Fortran code are provided in the file "instructions.txt". 

In addition to the code, we also provide a few example parameter files. These parameter files are used to define the survey parameters for the survey being used to detect a binary neutron star system of given orbital characteristics. The given parameter files were used in the merger rate analysis **Future prospects for ground-based gravitational wave detectors --- The Galactic double neutron star merger rate revisited** by **Pol, McLaughlin, and Lorimer** and can be used to calculate the degradation factors used in their dosurvey.py PsrPopPy source code (see the github repo for that paper [here](https://github.com/NihanPol/2018-DNS-merger-rate)). Note that Tobs, the integration time, will change based on the survey being used.

# Integration with Python

While the above code is written in Fortran, we have modified it so that it can be called into Python as a package and the degradation factor calculation can be called as a function in Python. If you do use this feature, please make sure to cite this package with the DOI given above.

This can be done using NumPy's f2py [functionality](https://docs.scipy.org/doc/numpy/f2py/). In this case, using the *.f90* files provided in this repository, do:

`python -m numpy.f2py -c xxx.f90 -m name_of_module`

where xxx.f90 is the code you want to work with (e.g. gamma1_real.f90) and "name_of_module" is the name of the python module corresponding to gamma1_real.f90. This will produce a "name_of_module.so" file.

Once that is done, in Python in the same folder where the "name_of_module.so" file lives, you can do
```
import numpy as np
import name_of_module #From the above command

#Make sure that the conversion has worked:
print(name_of_module.degfac.__doc__)
#This should print out a function signature and the input and output parameters

#This will now behave as a standard python function, where the 
different parameters can be passed through positional arguments.
#In case you miss a variable, the code will throw an error.

#As a test, try executing the function using the parameters 
given in one of the *.in files available in this repository.
```

If something does not work as advertised above, please get in touch with me and I will work on finding a fix.
