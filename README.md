# DNS-merger-rate-Update-2021

This is a repo containing all relevant files to peform the DNS merger rate estimation by [Pol et al. (2019)](https://arxiv.org/abs/1811.04086) with the extension including the lateral pulse profile of J1906+0746 as presented by [Grunthal et al. (2021)](https://arxiv.org/abs/2107.13307). 

# Using the files

To run the simulation as done in [Grunthal et al. (2021)](https://arxiv.org/abs/2107.13307), the following instruction might help
* Download and install the PsrPopPy simulation package, which can be found [here](https://github.com/devanshkv/PsrPopPy2)
* put the Python files contained in the *PsrPopPy_extension* folder into the lib/python folder in the PsrPopPy directory
* run the 'find_alphas.sh' script from the shell.
* Analyse the results using 'calculate_alpha_and_population_numbers.py' and 'calculate_merger_rates.py'. For the detailed use of the scripts, check out the original simulation and documentation under https://github.com/NihanPol/2018-DNS-merger-rate .

# Changes to the PsrPopPy files (depending only 1906+0746)

find_alpha: loads the profile data and creates spline
Survey: uses spline and creates rescaling factor for the beam amplitude

Make sure that all parameters required by the survey class are given in your simulation

