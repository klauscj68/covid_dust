This repository contains the codes used in our OSU campus effort to estimate
covid infection prevalence from environmental dust. It was written in Julia v1.6.0. 

### Installing the environment
The Julia environment needed to run the code is setup using the Project and
Manifest files. One can activate the environment by switching to the REPL
package manager and using the activate and instantiate commands. See for 
example
https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project.

Jupyter notebooks should automatically load installed packages without 
needing to activate the environment, although additional packages may need
to be installed. 

### Running MCMC
The julia files poimcmc.jl specifies most of the options, including
fixed parameter values, priors, number of samples, and how the
error is computed. Fixed parameter values are taken from data() in
poimcmc.jl. Priors are specified in mcmcrg() and logpi!() in poimcmc.jl.
mcmcrg() also specifies which parameters to vary.

To generate mcmc samples, one manually calls mcmcrun() with 
madatory argument specifying the desired the number of samples and optional
arguments specifying how often the chain is cycled before recording a sample
(for thinning), the random number generator state, the rate at which progress
through MCMC is printed to standard out, and whether the chain is continuing
an earlier MCMC run. At the conclusion of the routine, the script saves the
chain samples in MCMCsmp.csv and random number generator state in a series of
files of the form RNG*.csv. When restarting an mcmcrun(), these files should
be in the same folder as the directory being run from. Note that in that
case, mcmcrun() will overwrite those files, so you should manually save
backups if you want a larger chain.

### Postprocessing MCMC
The MCMC posterios are conveniently processed in the Jupyter notebook 
postprocess_mcmc.ipynb. It reads all the posterior samples from a single
MCMCsmp.csv file so multiple consecutive chain runs should be combined by the
user. 
