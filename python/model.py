#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is plotter for the SN Ia data by Ann Wang and Delilah Gates 
# Final Project for Physics 212: Cosmology 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


from __future__ import print_function #allow Python 2 to run print like Python 3
import sys, csv, random, math
from data_wrapper import data_stuff
from scipy import integrate as integral
import numpy as np 
import covariance as cov
import matplotlib.pyplot as plt
#import pylab as plt
import scipy.stats
from scipy.interpolate import spline


data_length= 0

c = 2.99792*np.power(10,5)
h = 0.7
H0 = 100 * h #km/s/Mpc

# fill in stuff here
alpha_mcmc =0.287259783104
beta_mcmc = 4.17179680554
Mp_mcmc = -19.360721486
dM_mcmc = -0.288884401992

om_mcmc =  0.468775219001

def dLumi(zhel, zcmb, om, ol, ok, func):
    # takes heliocentric redshift + CMB frame redshift values
    # and then outputs a luminosity distance!
    # arxiv 1601.01451 for reference, eq 2,3
    distance, distance_err = integral.quad(func, 0, zcmb, args=(om,ol,ok,))
    distance = (1+zhel)*c*distance/H0 
    return distance

def inverseHubble(z,om,ol,ok):
    # Calculates reduced hubble assuming lambda CDM
    # om is omega_matter
    return np.sqrt(om*np.power((1+z),3)+ok*np.power((1+z),2)+ol)

def M_calc(hostmass,Mp,dM):
    M = []
    for i in range(data_length):
        if hostmass[i] < 10:
            f = Mp
        else:
            f = Mp + dM
        M.append(f)
    return M

def mu_model_calc(zhel_data,zcmb_data,om,ol,ok):
    # calculates mu given z      
    dL = [dLumi(zhel_data[i],zcmb_data[i],om,ol,ok,inverseHubble) for i in range(data_length)]    
    mu_model=[25 + 5*np.log10(dL[i]) for i in range(data_length)]    
    return mu_model
                 
def mu_data_calc(m_B,x1,C,M,alpha,beta):
    # gives a functional form of mu given data + nuisance parameters       
    mu_data = [m_B[i]-M[i]+alpha*x1[i]-beta*C[i] for i in range(data_length)]    
    return mu_data    

def main():
    # call functions    
    data_wr = data_stuff()
    m_B_data, zcmb_data, zhel_data, c_data, x1_data, hostmass_data = data_wr.import_data("../data/jla_lcparams.txt")
    global data_length
    data_length = len(x1_data)
    print("Sample Size: ", data_length)
    M_data = M_calc(hostmass_data,Mp_mcmc,dM_mcmc)
    mu_data = mu_data_calc(m_B_data, x1_data, c_data, M_data, alpha_mcmc, beta_mcmc)
    mu_model = mu_model_calc(zhel_data,zcmb_data,om_mcmc,1-om_mcmc,0)
    error = cov.get_mu_cov(alpha_mcmc, beta_mcmc)
    mu_error = [np.sqrt(error[i][i]) for i in range(data_length)]
    z_error = [0 for i in range(data_length)]
    #    print(mu_error)
    print(len(error))
    print(len(mu_model))
    plt.rcParams.update({'font.size': 18})
#    mu_smooth = spline(np.asarray(zcmb_data),np.asarray(mu_model),np.asarray(z_smooth))
    plt.errorbar(zcmb_data,mu_data,yerr=mu_error,xerr=z_error,fmt="ro",markersize=5,label="SN Ia data")
    plt.plot(zcmb_data, mu_model, 'b.', linewidth=1,label="Theor. Pred.")
    legend = plt.legend(loc='lower right', shadow=False)
    plt.xlabel("Redshift")
    plt.ylabel(r"$\mu$")
    plt.axis([0,2,30.,48.])
    plt.savefig("../plots/mu.pdf")
    plt.clf()
if __name__ == "__main__":
    main()
