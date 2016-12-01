import sys, csv, random, math
from data_wrapper import data_stuff
from scipy import integrate as integral
import numpy as np
import covariance as cov

# final info from the mcmc

final_params_omega_l = []
final_params_omega_m = []
final_params_omega_k = []
final_l = []
final_weights = []

#global constants

c = 3*np.power(10,5)
h = 0.7
H0 = 100 * h #km/s/Mpc

# for now treating alpha and beta as constants that do not have to be fit for
# these are estimated values from Nielson et al.
alpha = 1
beta = 3.1
delta_m = -0.07
M_prime = -19

data_length = 0
nparam = 5
single_run_length=10
NUM_ITER = 10


def mc_standard_err():
    # function for calculating mcmc error
    # TBD
    return 0

def om_error():
    # function for calculating om error
    # TBD
    return 0

def likelihood_calc(mu_model,mu_data):
    # evaluates likelihood given the model
    # basic error
    likelihood=[]
    P = 0    
    n = range(data_length)
    for i in n:
        p = (mu_data[i]-mu_model[i])*(mu_data[i]-mu_model[i])/mu_model[i]
        print "DIFF: ",mu_data[i]-mu_model[i]
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
#    om = random.uniform(0.2,0.4)
    ol = 1 - om
    ok = 0
    #ok=floor(rand*3)-1
    return om, ol, ok

def nuisance_draw():
    global alpha
    global beta
    global delta_m
    global M_prime
    alpha=random.uniform(0,1)
    alpha=random.gauss(0.141,0.006)
    beta=random.uniform(2.5,5)
    beta = random.gauss(3.101,0.075)
    delta_m=random.uniform(-1,0)
    delta_m = random.gauss(-0.070,0.023)
    M_prime = random.uniform(-20,-18)
    M_prime = random.gauss(-19.05,0.02)
    # for gaussian random.gauss(mu, sigma)
    return alpha, beta, delta_m, M_prime

def dLumi(zhel, zcmb, om, ol, ok, func):
    # takes heliocentric redshift + CMB frame redshift values
    # and then outputs a luminosity distance!
    # arxiv 1601.01451 for reference, eq 2,3
    distance, distance_err = integral.quad(func, 0, zcmb, args=(om,ol,ok,))
    distance = (1+zhel)*c*distance/H0
#    print "distance: ", distance, " zhel: ", zhel
    return distance

def inverseHubble(z,om,ol,ok):
    # Calculates reduced hubble assuming lambda CDM
    # om is omega_matter
    return np.sqrt(om*np.power((1+z),3)+ok*np.power((1+z),2)+ol)

def M_step(host_mass,delta_m,M_prime):
    M = []
    for i in range(data_length):
        if host_mass[i] < 10:
            f = M_prime
        else:
            f = M_prime + delta_m
        M.append(f)
    return M

def mu_model_calc(zhel_data,zcmb_data,omega_m,omega_l,omega_k):
    # calculates mu given z
    mu_model=[];

    for i in range(data_length):
        dL = dLumi(zhel_data[i],zcmb_data[i],omega_m,omega_l,omega_k,inverseHubble)
        mu_temp = 25+5*np.log10(dL)
        mu_model.append(mu_temp)
#        print "model: ",mu_temp
    return mu_model

def mu_data_calc(m_B,x1,C,M_list,alpha,beta):
    # gives a functional form of mu given data + nuisance parameters
    mu_data=[]
    
    n = range(data_length)
    for i in n:
        f=m_B[i]-M_list[i]+alpha*x1[i]-beta*C[i]
        mu_data.append(f)
#        print "data: " ,f
    return mu_data

def likelihood_draw(zhel_data,zcmb_data,m_B_data,x1_data,c_data,host_data):
    om,ol,ok = omega_draw()
    dM = delta_m
    a = alpha
    b = beta
    Mp = M_prime
#    a,b,dM,Mp = nuisance_draw()
    M_list = M_step(host_data,dM,M_prime)
    #print "calculating model mu"
    mu = mu_model_calc(zhel_data,zcmb_data,om,ok,ol)
    #print "calculating estimated mu"
    mu_hat = mu_data_calc(m_B_data,x1_data,c_data,M_list,alpha,beta)
    P = likelihood_calc(mu_hat, mu)
    return P, mu, mu_hat, om, ol, ok, a, b, dM, Mp

def step_mcmc(single_run_length,zhel_data,zcmb_data,m_B_data,x1_data,c_data,host_data):
    # step through the mcmc
    P = 0
    mu = 0
    mu_hat = 0
    om = 0
    ol = 0
    ok = 0
    a = 0
    b = 0
    dM = 0
    Mp = 0
    w = 0 # weight
    
    P, mu, mu_hat, om, ol, ok, a, b, dM, Mp = likelihood_draw(zhel_data,zcmb_data,m_B_data,x1_data,c_data,host_data);
    #for i in range(single_run_length):
    while(True):
        print "OLD OM: ", om
        P_new, mu_temp, mu_hat_temp, om_temp, ol_temp, ok_temp, a_temp, b_temp, dM_temp, Mp_temp = likelihood_draw(zhel_data,zcmb_data,m_B_data,x1_data,c_data,host_data)
        print "NEW OM: ", om_temp
        r = P_new/P
        print "ratio: ", r
        if r > 1:
            w = w+1
        else:
            w = 0
            P = P_new
            om = om_temp
            ol = ol_temp
            ok = ok_temp
            a = a_temp
            b = b_temp
            dM = dM_temp
            Mp = Mp_temp
        # if (P/(data_length-nparam) < 2):
        #      print "ooo, small x2!"
        #      break
        print "dof: ",data_length-nparam
        print "Likelihood: ",P/(data_length-nparam)
        print "alpha, beta, dM, Mp: ",a,b,dM,Mp
    return om, ol, ok, w, P
    
def covariance():
    # think about this later
    return 0

def main():
    # call functions
    data_wr = data_stuff()
    m_B_data, zcmb_data, zhel_data, c_data, x1_data, host_data = data_wr.import_data("../data/jla_lcparams.txt")
    global data_length
    data_length = len(x1_data)
    print "Sample Size: ", data_length
    n = range(NUM_ITER)
    print cov.get_mu_cov(0.13, 3.1)
    for i in n:
        print "Iteration: ", i+1
        param_omega_m, param_omega_l, param_omega_k, weight, likelihood = step_mcmc(single_run_length,zhel_data,zcmb_data,m_B_data,x1_data,c_data,host_data)
        final_params_omega_m.append(param_omega_m)
        final_params_omega_l.append(param_omega_l)
        final_params_omega_k.append(param_omega_k)
        final_l.append(likelihood)
        final_weights.append(weight)
    print final_params_omega_m,final_params_omega_l,final_params_omega_k, final_l, final_weights
    weighted_omega_m_avg = 0
    total_weight = 0
    for i in range(len(final_params_omega_m)):
        weighted_omega_m_avg = weighted_omega_m_avg + final_params_omega_m[i]*final_weights[i]
        total_weight = total_weight + final_weights[i]
    print "WEIGHTED OMEGA_M AVG: ", weighted_omega_m_avg/total_weight
    print "alpha, beta, M', dM", alpha, beta, M_prime, delta_m
if __name__ == "__main__":
    main()
