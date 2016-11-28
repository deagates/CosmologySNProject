import sys, csv, random, math
from data_wrapper import data_stuff
from scipy import integrate as integral
import numpy as np 

# final info from the mcmc

final_params_omega_l = []
final_params_omega_m = []
final_params_omega_k = []
final_l = []
final_weights = []

#global constants

c = 1
h = 0.673
H0 = 100 * h #km/s/Mpc

# for now treating alpha and beta as constants that do not have to be fit for
# these are estimated values from Nielson et al.
alpha = 0.1
beta = 3
M = -19

data_length = 0

single_run_length=10
NUM_ITER = 10

def likelihood_calc(mu_model,mu_data):
    # evaluates likelihood given the model
    # basic error
    likelihood=[]
    P = 0    
    n = range(data_length-1)
    for i in n:
        p = (mu_data[i]-mu_model[i])*(mu_data[i]-mu_model[i])/mu_model[i]
        P = P + p;
    return P


def omega_draw():
    om = random.random()
    ol = 1 - om
    ok = 1
    #ok=floor(rand*3)-1
    return om, ol, ok

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
    
    def ii(x):
        return H0/H(x,omega_m,omega_l,omega_k)
    
    n = range(data_length-1)
    for i in n:
        #z=M[i]
        #I, err=sc.integrate.quad(ii(x),0,z)
        #dL=(1+z)*(dH/sqrt(omega_k))*math.asin(sqrt(omgea_k)*I)
        dL = dLumi(zhel_data[i],zcmb_data[i],omega_m,omega_l,omega_k,inverseHubble)
        mu_temp = 25 + 5*np.log10(dL)
        mu_model.append(mu_temp)
    return mu_model

def mu_data_calc(m_B,x1,C):
    # gives a functional form of mu given data + nuisance parameters
    mu_data=[]
    
    n = range(data_length-1)
    for i in n:
        f=m_B[i]-M+alpha*x1[i]-beta*C[i]
        mu_data.append(f)
    return mu_data


def mu_error():
    return 0
    # gives error associated
    

def likelihood_draw(zhel_data,zcmb_data,m_B_data,x1_data,c_data):
    om,ol,ok = omega_draw()
    #print "calculating model mu"
    mu = mu_model_calc(zhel_data,zcmb_data,om,ok,ol)
    #print "calculating estimated mu"
    mu_hat = mu_data_calc(m_B_data,x1_data,c_data)
    P = likelihood_calc(mu_hat, mu)
    return P, mu, mu_hat, om, ol, ok


def step_mcmc(single_run_length,zhel_data,zcmb_data,m_B_data,x1_data,c_data):
    # evaluates one step of mcmc algorithm
    # put weight here for now
    #initialize 
    P_old, mu, mu_hat, om, ol, ok = likelihood_draw(zhel_data,zcmb_data,m_B_data,x1_data,c_data);
#    print P_old, mu, mu_hat, om, ol, ok
    w=0
    n = range(single_run_length)
    for i in n:
        P_new, mu, mu_hat, om, ol, ok = likelihood_draw(zhel_data,zcmb_data,m_B_data,x1_data,c_data)
        r=P_new/P_old;
        if r>1:
            w=w+1;
        else:
            w=0;
            P_old=P_new;
    return om, ol, ok, w, P_new

    
def covariance():
    # think about this later
    return 0

def main():
    # call functions
    data_wr = data_stuff()
    m_B_data, zcmb_data, zhel_data, c_data, x1_data = data_wr.import_data("../data/jla_lcparams.txt")
    global data_length
    data_length = len(x1_data)
    print "Sample Size: ", data_length
    n = range(NUM_ITER)
    for i in n:
        print "Iteration: ", i+1
        param_omega_m, param_omega_l, param_omega_k, weight, likelihood = step_mcmc(single_run_length,zhel_data,zcmb_data,m_B_data,x1_data,c_data)
        final_params_omega_m.append(param_omega_m)
        final_params_omega_l.append(param_omega_l)
        final_params_omega_k.append(param_omega_k)
        final_l.append(likelihood)
        final_weights.append(weight)
    #print final_params_omega_m,final_params_omega_l,final_params_omega_k, final_l, final_weights
    weighted_omega_m_avg = 0
    total_weight = 0
    for i in range(len(final_params_omega_m)):
        weighted_omega_m_avg = weighted_omega_m_avg + final_params_omega_m[i]*final_weights[i]
        total_weight = total_weight + final_weights[i]
    print "WEIGHTED OMEGA_M AVG: ", weighted_omega_m_avg/total_weight

if __name__ == "__main__":
    main()
