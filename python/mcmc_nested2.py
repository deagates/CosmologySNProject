from __future__ import print_function #allow Python 2 to run print like Python 3
import sys, csv, random, math
from data_wrapper import data_stuff
from scipy import integrate as integral
import numpy as np 
import covariance as cov
#import matplotlib.pyplot as plt
import pylab as plt

# final info from the mcmc

final_params_omega_l = []
final_params_omega_m = []
final_params_omega_k = []

final_params_alpha = []
final_params_beta = []
final_params_Mp = []
final_params_dM = []

final_l = []
final_weights = []

#global constants

c = 3*np.power(10,5)
h = 0.673
H0 = 100 * h #km/s/Mpc
plot_mcmc_mini_conv = "off"
plot_mcmc_conv = "off"
user = "D"

# for now treating alpha and beta as constants that do not have to be fit for
# these are estimated values from Nielson et al.
alpha = 0.1
beta = 3
#M_prime = -19.03

data_length = 0

single_run_length=10
single_run_length_mini=15
NUM_ITER = 1

def likelihood_calc(mu_model,mu_data):
    # evaluates likelihood given the model
    # basic error
    P = 0    
    n = range(data_length)
    for i in n:
        p = (mu_data[i]-mu_model[i])*(mu_data[i]-mu_model[i])/mu_model[i]
        P = P + p;
    return P

def X2_likelihood_calc(mu_model,mu_data):
    # evaluates X2 likelihood given the model
    mu_data_np = np.array(mu_data)
#    print np.shape(mu_data_np)
    mu_model_np = np.array(mu_model)
#    print np.shape(mu_model_np)
    mu_diff = mu_data_np - mu_model_np
    inv_covmat = np.linalg.inv(cov.get_mu_cov(alpha, beta))
    return mu_diff.dot(inv_covmat.dot(mu_diff))

def omega_draw():
    om = random.random()
    ol = 1 - om
    ok = 0
    #ok=floor(rand*3)-1
    return om, ol, ok

def nuisance_draw():
    global alpha
    global beta
    global delta_m
    global M_prime
    alpha=random.uniform(0,4)
    #alpha=random.gauss(0.141,0.006)
    beta=random.uniform(1,5)
#    beta = random.gauss(3.101,0.075)
    delta_m=random.uniform(-1,0)
#    delta_m = random.gauss(-0.070,0.023)
    M_prime = random.uniform(-20,-18)
#    M_prime = random.gauss(-19.05,0.02)
    # for gaussian random.gauss(mu, sigma)
    return alpha, beta, delta_m, M_prime
    
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

def mu_model_calc(zhel_data,zcmb_data,omega_m,omega_l,omega_k):
    # calculates mu given z
    mu_model=[];
    
    #def ii(x):
    #    return H0/H(x,omega_m,omega_l,omega_k)
    
    n = range(data_length)
    for i in n:
        #z=M[i]
        #I, err=sc.integrate.quad(ii(x),0,z)
        #dL=(1+z)*(dH/sqrt(omega_k))*math.asin(sqrt(omgea_k)*I)
        dL = dLumi(zhel_data[i],zcmb_data[i],omega_m,omega_l,omega_k,inverseHubble)
        mu_temp = 25 + 5*np.log10(dL)
        mu_model.append(mu_temp)
    return mu_model

def M_calc(hostmass,Mp,dM):
    M = []
    n = range(data_length)
    for i in n:
        if hostmass[i] < np.power(10,10):
            f=Mp
        else:
            f=Mp+dM
        M.append(f)
    return M
                 
def mu_data_calc(m_B,x1,C,M,alpha,beta):
    # gives a functional form of mu given data + nuisance parameters
    mu_data=[]
    
    n = range(data_length)
    for i in n:        
        f=m_B[i]-M[i]+alpha*x1[i]-beta*C[i]
        mu_data.append(f)
    return mu_data


def mu_error():
    return 0
    # gives error associated
    
    
def likelihood_draw_mini(m_B_data,x1_data,c_data,hostmass,mu):    
    alpha, beta, Mp, dM = nuisance_draw()
    M=M_calc(hostmass,Mp,dM)    
    mu_hat = mu_data_calc(m_B_data,x1_data,c_data,M,alpha,beta)
    P = X2_likelihood_calc(mu_hat, mu)
    return P, mu_hat, alpha, beta, Mp, dM
    
    
def step_mcmc_mini(m_B_data,x1_data,c_data,hostmass,mu) : 
    P_mini, mu_hat, alpha, beta, Mp, dM = likelihood_draw_mini(m_B_data,x1_data,c_data,hostmass,mu)
    #print("pmini:",P_mini)
    
    if plot_mcmc_mini_conv == "on":
        Ppmini= []
        Ppmini_temp = []
        Rr = []
        Ii = []
        
    
    w_mini=0
    n = range(single_run_length_mini)
    
    for i in n:
        P_mini_temp, mu_hat_temp, alpha_temp, beta_temp , Mp_temp, dM_temp = likelihood_draw_mini(m_B_data,x1_data,c_data,hostmass,mu)
        r=P_mini_temp/P_mini
          
        if plot_mcmc_mini_conv == "on":
            Ii.append(i)
            Ppmini.append(P_mini/10000+.02)
            Ppmini_temp.append(P_mini_temp/10000+.02)
            Rr.append(r)
 
        if r>1:
            w_mini=w_mini+1;
        else:
            w_mini=0;
            P_mini=P_mini_temp;
            mu_hat=mu_hat_temp;
            alpha=alpha_temp;
            beta=beta_temp;
            Mp=Mp_temp
            dM=dM_temp;
            
    if plot_mcmc_mini_conv == "on":   
        Aa=[1]* single_run_length_mini    
        plt.figure(1)
        plt.subplot(211)
        plt.plot(Ii,Ppmini,'b',label='P')
        plt.plot(Ii,Ppmini_temp,'g+',label='P_temp')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.ylim((0,Ppmini[0]+10))
        plt.ylabel('X2_likelihood')
        plt.subplot(212)
        plt.plot(Ii,Rr,'r*', label='r from mcmc')
        plt.plot(Ii,Aa,'r', label = 'r=1')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.ylim((0,10))
        plt.ylabel('r=P_temp/P')
        plt.xlabel('mcmc step')
        plt.show()
    return alpha, beta, Mp,  dM, w_mini, P_mini, mu_hat    

def likelihood_draw(zhel_data,zcmb_data,m_B_data,x1_data,c_data,hostmass):
    om,ol,ok = omega_draw()
    mu = mu_model_calc(zhel_data,zcmb_data,om,ok,ol)
    alpha, beta, Mp, dM, w_mini, P_mini, mu_hat=step_mcmc_mini(m_B_data,x1_data,c_data,hostmass,mu) 
    P = X2_likelihood_calc(mu_hat, mu)
    return P, mu, mu_hat, om, ol, ok, alpha, beta, Mp, dM
    
def step_mcmc(single_run_length,zhel_data,zcmb_data,m_B_data,x1_data,c_data,hostmass):
    # evaluates one step of mcmc algorithm
    # put weight here for now
    #initialize 
    P, mu, mu_hat, om, ol, ok,alpha, beta, Mp, dM = likelihood_draw(zhel_data,zcmb_data,m_B_data,x1_data,c_data,hostmass)
#   print P_old, mu, mu_hat, om, ol, ok

    if plot_mcmc_conv == "on":
        Pp= []
        Pp_temp = []
        Rr = []
        Ii = []


    w=0
    n = range(single_run_length)
    for i in n:
        P_temp, mu_temp, mu_hat_temp, om_temp, ol_temp, ok_temp ,alpha_temp, beta_temp, Mp_temp, dM_temp = likelihood_draw(zhel_data,zcmb_data,m_B_data,x1_data,c_data,hostmass)
        r=P_temp/P;
        
        if plot_mcmc_conv == "on":
            Ii.append(i)
            Pp.append(P/10000-0.1)
            Pp_temp.append(P_temp/10000-0.1)
            Rr.append(r)
        
        if r>1:
            w=w+1
        else:
            w=0;
            P=P_temp
            om=om_temp
            ok=ok_temp
            om=om_temp
            #mu_hat=mu_hat_temp
            alpha=alpha_temp
            beta=beta_temp
            Mp=Mp_temp
            dM=dM_temp
        
    if plot_mcmc_conv == "on":   
        Aa=[1]* single_run_length 
        #plt.figure(1)
        #plt.subplot(211)
        #plt.plot(Ii,Pp,'b',Ii,Pp_temp,'g+')
        #plt.ylim((0,Pp[0]+10))
        #plt.subplot(212)
        #plt.plot(Ii,Rr,'r*',Ii,Aa,'r')
        #plt.ylim((0,10))
        plt.figure(1)
        plt.subplot(211)
        plt.plot(Ii,Pp,'b',label='P')
        plt.plot(Ii,Pp_temp,'g+',label='P_temp')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.ylim((0,Pp[0]+10))
        plt.ylabel('X2_likelihood')
        plt.subplot(212)
        plt.plot(Ii,Rr,'r*', label='r from mcmc')
        plt.plot(Ii,Aa,'r', label = 'r=1')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.ylim((0,10))
        plt.ylabel('r=P_temp/P')
        plt.xlabel('mcmc step')
        plt.show()  
        
    return om, ol, ok,alpha, beta, Mp, dM, w, P

    
def covariance():
    # think about this later
    return 0

def main():
    # call functions    
    data_wr = data_stuff()
    m_B_data, zcmb_data, zhel_data, c_data, x1_data, hostmass_data = data_wr.import_data("../data/jla_lcparams.txt")
    global data_length
    data_length = len(x1_data)
    print("Sample Size: ", data_length)
    n = range(NUM_ITER)
    #print(cov.get_mu_cov(0.13, 3.1))
    for i in n:
        print("Iteration: ", i+1)
        param_omega_m, param_omega_l, param_omega_k, alpha, beta, Mp, dM, weight, likelihood = step_mcmc(single_run_length,zhel_data,zcmb_data,m_B_data,x1_data,c_data,hostmass_data)
        final_params_omega_m.append(param_omega_m)
        final_params_omega_l.append(param_omega_l)
        final_params_omega_k.append(param_omega_k)
        final_params_alpha.append(alpha)
        final_params_beta.append(beta)
        final_params_Mp.append(Mp)
        final_params_dM.append(dM)
        final_l.append(likelihood)
        final_weights.append(weight)
        print("~~~")       
    print(final_params_omega_m,final_params_omega_l,final_params_omega_k, final_l, final_weights) 
    weighted_omega_m_avg = 0
    total_weight = 0
    for i in range(len(final_params_omega_m)):
        weighted_omega_m_avg = weighted_omega_m_avg + final_params_omega_m[i]*final_weights[i]
        total_weight = total_weight + final_weights[i]
    print("WEIGHTED OMEGA_M AVG: ", weighted_omega_m_avg/total_weight)
    print("alpha, beta, M', dM", alpha, beta, M_prime, delta_m)

if __name__ == "__main__":
    main()
