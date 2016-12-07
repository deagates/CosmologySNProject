
# Likelihood analysis of SN Ia Data #

Everything in the data folder is taken from Betoule et al., found on http://supernovae.in2p3.fr/sdss_snls_jla/ReadMe.html. This includes the JLA catalogue of SN Ia, along with the covariance matrix files and python converter (given \alpha, \beta) found in /data/covmat/.

The main scripts for the analysis are in /python/, where Run_mcmc_extenral.py calls the main file mcmc.py. The data is inputed using data_wrapper.py, and covariance.py contains the aforementioned python covariance matrix converter provided by Betoule et al. 
