# To be run by "Run_mcmc_external.py"
# So multiple runs can be preformed changing mcmc run lengths to check convergence
#and so data from multiple runs can easily be saved into different files
#mcmc_nested_external.py cannot be run of its own accord

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#This is the MCMC progam by Ann Wang and Delilah Gates 
#final project for Phys212: Cosmomology 
#We take in SNe Ia data and run an MCMC to find the most likely values for 
#Omega_m and Omega_Lambda, the matter and dark energy content of the universe
#This was inspired by a recent Nature paper (http://www.nature.com/articles/srep35596) 
# which claims SNe data only suggests maginally that there is dark energy 
#in the universe
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


from __future__ import print_function #allow Python 2 to run print like Python 3
import sys, csv, random, math
from data_wrapper import data_stuff
from scipy import integrate as integral
import numpy as np 
import covariance as cov
#import matplotlib.pyplot as plt
import pylab as plt



data_length= 0
nparam = 5 # number of parameters that can vary (omega_m, and 4 nusicene parameters)
#single_run_length=3
#single_run_length_mini=3
#NUM_ITER = 3
#ResultsName='result.txt'

#global constants
c = 3*np.power(10,5)
h = 0.7
H0 = 100 * h #km/s/Mpc
plot_mcmc_mini_conv = "off"
plot_mcmc_conv = "off"
user = "D"


# intermediate iterations will be saved here
alpha_save = ['alpha']
beta_save = ['beta']
Mp_save = ['Mp']
dM_save = ['dM']
ol_save = ['omega_Lambda']
om_save = ['omega_m']
l_save = ['X2_likelihood']
w_save = ['weight']

#final parameters from each full mcmc iteration saved here
final_params_omega_l = []
final_params_omega_m = []
final_params_omega_k = []

final_params_alpha = []
final_params_beta = []
final_params_Mp = []
final_params_dM = []

final_l = []
final_weights = []


#==============================================================================
#obsolete!: naive likelihood parameter before we included errors
# def likelihood_calc(mu_model,mu_data):
#     # evaluates likelihood given the model
#     # basic error
#     P = 0    
#     for i in range(data_length):
#         p = (mu_data[i]-mu_model[i])*(mu_data[i]-mu_model[i])/mu_model[i]
#         P = P + p;
#     return P
#==============================================================================

def X2_likelihood_calc(mu_model,mu_data):
    # evaluates X2 likelihood given the model
    mu_data_np = np.array(mu_data)
#    print np.shape(mu_data_np)
    mu_model_np = np.array(mu_model)
#    print np.shape(mu_model_np)
    mu_diff = mu_data_np - mu_model_np
    inv_covmat = np.linalg.inv(cov.get_mu_cov(alpha, beta))
    return mu_diff.dot(inv_covmat.dot(mu_diff))/(data_length-nparam)

def omega_draw():
    om = random.random() #omega_m (matter content is between 0 and 1, flat prior)
    ol = 1 - om  #omega_lambda (dark energy content)
    ok = 0 #for simplicity we impose omega_k= 0 (flat universe)
    #ok=floor(rand*3)-1
    return om, ol, ok

def nuisance_draw():
    # this function randomly draws nuisance paramteres
    # these are the 4 nuisance parameters the SNe luminosoty magnitude model will
    # use to come up with the expected value of mu_hat based on the red shift

    global alpha
    global beta
    global delta_m
    global M_prime

    alpha=random.uniform(0,4)
    #alpha=random.gauss(0.141,0.006)
    beta=random.uniform(1,5)
    #beta = random.gauss(3.101,0.075)
    delta_m=random.uniform(-1,0)
    #delta_m = random.gauss(-0.070,0.023)
    M_prime = random.uniform(-20,-18)
    #M_prime = random.gauss(-19.05,0.02)
    # for gaussian random.gauss(mu, sigma)    

    # below are estimated values from Nielson et al.
    #alpha = 0.1
    #beta = 3
    #M_prime = -19.03
    return alpha, beta, M_prime, delta_m
    
    
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
            f = M_prime
        else:
            f = M_prime + delta_m
        M.append(f)
        
    return M

def mu_model_calc(zhel_data,zcmb_data,om,ol,ok):
    # calculates mu given z
#==============================================================================
#     mu_model=[];
#     
#     for i in range(data_length):
#         dL = dLumi(zhel_data[i],zcmb_data[i],omega_m,omega_l,omega_k,inverseHubble)
#         mu_temp = 25 + 5*np.log10(dL)
#         mu_model.append(mu_temp)
#==============================================================================
        
    dL = [dLumi(zhel_data[i],zcmb_data[i],om,ol,ok,inverseHubble) for i in range(data_length)]    
    mu_model=[25 + 5*np.log10(dL[i]) for i in range(data_length)]    
    return mu_model
                 
def mu_data_calc(m_B,x1,C,M,alpha,beta):
    # gives a functional form of mu given data + nuisance parameters
#==============================================================================
#     mu_data=[]
#     
#     for i in range(data_length):        
#         f=m_B[i]-M[i]+alpha*x1[i]-beta*C[i]
#         mu_data.append(f)
#==============================================================================
        
    mu_data=[m_B[i]-M[i]+alpha*x1[i]-beta*C[i] for i in range(data_length)]    
    return mu_data
    
    
def likelihood_draw_mini(m_B_data,x1_data,c_data,hostmass,mu):
    # draws nuisance parameters, calculated mu_hat. compares a given mu  to 
    # mu_hat and calulates X2 likelihood  
    
    #draw nusiance pareters
    alpha, beta, Mp, dM = nuisance_draw() 
    #print(alpha)
    
    #calculate mu_hat
    M = M_calc(hostmass,Mp,dM)    
    mu_hat = mu_data_calc(m_B_data,x1_data,c_data,M,alpha,beta)
    
    #calculate X2_likelihood
    P = X2_likelihood_calc(mu_hat, mu)
    return P, mu_hat, alpha, beta, Mp, dM    
    
def step_mcmc_mini(single_run_length_mini,m_B_data,x1_data,c_data,hostmass,mu) : 
    # runs mcmc over nuisance parameters for some given mu,to converge on smaller
    # X2 likelihood
    
    #initalize
    #draws nusiance parameters, and calulated mu_hat, and X2_likelihood for given mu
    #base step: with likelihood and corresponding parameters stored 
    P_mini, mu_hat, alpha, beta, Mp, dM = likelihood_draw_mini(m_B_data,x1_data,c_data,hostmass,mu)
    #print("pmini:",P_mini)
#==============================================================================
#     if plot_mcmc_mini_conv == "on":
#         Ppmini = []
#         Ppmini_temp = []
#         Rr = []
#         Ii = []
#==============================================================================
    #initalize weight    
    w_mini = 0
    
    for i in range(single_run_length_mini):
#==============================================================================
#        if np.mod(i,100) == 0:
#            print("iterating. Mini Draws left",single_run_length_mini-i)
#         else:
#             print("iterating...")
#==============================================================================
        
        #draws nusiance parameters, and calulated mu_hat, and X2_likelihood for given mu
        #temporary     
        P_mini_temp, mu_hat_temp, alpha_temp, beta_temp, Mp_temp, dM_temp = likelihood_draw_mini(m_B_data,x1_data,c_data,hostmass,mu)
        #compares pervious stored X2 likelihood to new temporary X2 likelhood
        r = P_mini_temp/P_mini
          
#==============================================================================
#         if plot_mcmc_mini_conv == "on":
#             Ii.append(i)
#             Ppmini.append(P_mini)
#             Ppmini_temp.append(P_mini_temp)
#             Rr.append(r)
#==============================================================================
        # if new temporary X2 likelhood is > than pervious stored X2 likelihood... 
        if r > 1: 
            # weight goe up by one
            w_mini = w_mini + 1; 
        # otherwise (new temporary X2 likelhood is < or =  pervious stored X2 likelihood...)   
        else:
            #set weight back to zero
            w_mini = 0;
            #and set X2 likelihood and corresponding parameters stored for future comparison
            P_mini = P_mini_temp;
            mu_hat = mu_hat_temp;
            alpha = alpha_temp;
            beta = beta_temp;
            Mp = Mp_temp
            dM = dM_temp;
            
#==============================================================================
#     Code for ploting convergence to make sure procees is working correctly  
#     if plot_mcmc_mini_conv == "on":   
#         Aa = [1] * single_run_length_mini    
#         plt.figure(1)
#         plt.subplot(211)
#         plt.plot(Ii,Ppmini,'b',label='P')
#         plt.plot(Ii,Ppmini_temp,'g+',label='P_temp')
#         plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#         plt.ylim((0,Ppmini[0]+10))
#         plt.ylabel('X2_likelihood')
#         plt.subplot(212)
#         plt.plot(Ii,Rr,'r*', label='r from mcmc')
#         plt.plot(Ii,Aa,'r', label = 'r = 1')
#         plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#         plt.ylim((0,10))
#         plt.ylabel('r=P_temp/P')
#         plt.xlabel('mcmc step')
#         #plt.show()
#==============================================================================
    return alpha, beta, Mp,  dM, w_mini, P_mini, mu_hat    

def likelihood_draw(single_run_length_mini,zhel_data,zcmb_data,m_B_data,x1_data,c_data,hostmass):
    # draws parameters, calculates mu_hat and mu. and compares them to calulates X2 likelihood  
    
    #draws omegas
    om,ol,ok = omega_draw()
    #calculates mu
    mu = mu_model_calc(zhel_data,zcmb_data,om,ok,ol)
    print("drew different omegas")
    # draws nusiance parameters and calculates mu_hat
    alpha, beta, Mp, dM, w_mini, P_mini, mu_hat = step_mcmc_mini(single_run_length_mini,m_B_data,x1_data,c_data,hostmass,mu) 
    # calculates X2 likelihood from mu and mu_hat
    P = X2_likelihood_calc(mu_hat, mu)
    return P, mu, mu_hat, om, ol, ok, alpha, beta, Mp, dM
    
def step_mcmc(single_run_length,single_run_length_mini,zhel_data,zcmb_data,m_B_data,x1_data,c_data,hostmass):
    # evaluates one step of mcmc algorithm
    
    
    #initialize 
    #draws parameters, and calulates mu and mu_hat, calculates X2_likelihood 
    #base step: with likelihood and corresponding parameters stored 
    P, mu, mu_hat, om, ol, ok, alpha, beta, Mp, dM = likelihood_draw(single_run_length_mini,zhel_data,zcmb_data,m_B_data,x1_data,c_data,hostmass)
#   print P_old, mu, mu_hat, om, ol, ok

    #All data is going to be stored in these global functions
    global alpha_save
    global beta_save
    global dM_save
    global Mp_save
    global om_save
    global ol_save
    global l_save
    global w_save
    
    alpha_save.append('new iteration')
    beta_save.append('~')
    dM_save.append('~')
    Mp_save.append('~')
    om_save.append('~')
    ol_save.append('~')
    l_save.append('~')
    w_save.append('~')

    if plot_mcmc_conv == "on":
        Pp= []
        Pp_temp = []
        Rr = []
        Ii = []

    w=0
    for i in range(single_run_length):        
        alpha_save.append(alpha)
        beta_save.append(beta)
        dM_save.append(dM)
        Mp_save.append(Mp)
        om_save.append(om)
        ol_save.append(ol)
        l_save.append(P)
        w_save.append(w)
        
        P_temp, mu_temp, mu_hat_temp, om_temp, ol_temp, ok_temp ,alpha_temp, beta_temp, Mp_temp, dM_temp = likelihood_draw(single_run_length_mini,zhel_data,zcmb_data,m_B_data,x1_data,c_data,hostmass)
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
            ol=ol_temp
            #mu_hat=mu_hat_temp
            alpha=alpha_temp
            beta=beta_temp
            Mp=Mp_temp
            dM=dM_temp
        
#==============================================================================
#     if plot_mcmc_conv == "on":   
#         Aa=[1]* single_run_length 
#         plt.figure(1)
#         plt.subplot(211)
#         plt.plot(Ii,Pp,'b',label='P')
#         plt.plot(Ii,Pp_temp,'g+',label='P_temp')
#         plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#         plt.ylim((0,Pp[0]+10))
#         plt.ylabel('X2_likelihood')
#         plt.subplot(212)
#         plt.plot(Ii,Rr,'r*', label='r from mcmc')
#         plt.plot(Ii,Aa,'r', label = 'r=1')
#         plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#         plt.ylim((0,10))
#         plt.ylabel('r=P_temp/P')
#         plt.xlabel('mcmc step')
#         #plt.show()  
#         plt.savefig("mini.pdf")
#==============================================================================
        
    return om, ol, ok,alpha, beta, Mp, dM, w, P

def main(NUM_ITER,single_run_length,single_run_length_mini,ResultsName):
    # call functions    
    data_wr = data_stuff()
    m_B_data, zcmb_data, zhel_data, c_data, x1_data, hostmass_data = data_wr.import_data("../data/jla_lcparams.txt")
    global data_length
    data_length = len(x1_data)
    print("Sample Size: ", data_length)
    #print(cov.get_mu_cov(0.13, 3.1))
    for i in range(NUM_ITER):
        print("Iterations left: ", NUM_ITER-i)
        param_omega_m, param_omega_l, param_omega_k, alpha, beta, Mp, dM, weight, likelihood = step_mcmc(single_run_length,single_run_length_mini,zhel_data,zcmb_data,m_B_data,x1_data,c_data,hostmass_data)
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
    #print(final_params_omega_m,final_params_omega_l,final_params_omega_k, final_l, final_weights) 
#==============================================================================
#     weighted_omega_m_avg = 0
#     total_weight = 0
#     for i in range(len(final_params_omega_m)):
#         weighted_omega_m_avg = weighted_omega_m_avg + final_params_omega_m[i]*final_weights[i]
#         total_weight = total_weight + final_weights[i]
#==============================================================================
    
    weighted_omega_m_avg=sum([final_params_omega_m[i]*final_weights[i] for i in range(len(final_params_omega_m))])    
    weighted_omega_l_avg=sum([final_params_omega_l[i]*final_weights[i] for i in range(len(final_params_omega_m))])
    weighted_alpha_avg=sum([final_params_alpha[i]*final_weights[i] for i in range(len(final_params_omega_m))])
    weighted_beta_avg=sum([final_params_beta[i]*final_weights[i] for i in range(len(final_params_omega_m))])
    weighted_Mp_avg=sum([final_params_Mp[i]*final_weights[i] for i in range(len(final_params_omega_m))])
    weighted_dM_avg=sum([final_params_dM[i]*final_weights[i] for i in range(len(final_params_omega_m))])
    
    total_weight=sum([final_weights[i] for i in range(len(final_params_omega_m))])    
    
    print("WEIGHTED OMEGA_M AVG: ", weighted_omega_m_avg/total_weight)
    print("WEIGHTED OMEGA_l AVG: ", weighted_omega_l_avg/total_weight)
    print("WEIGHTED Alpha AVG: ", weighted_alpha_avg/total_weight)
    print("WEIGHTED Beta AVG: ", weighted_beta_avg/total_weight)
    print("WEIGHTED Mp AVG: ", weighted_Mp_avg/total_weight)
    print("WEIGHTED dM AVG: ", weighted_dM_avg/total_weight)
    
    print("omega_, omega_l:",param_omega_m, param_omega_l)
    print("alpha, beta:", alpha, beta)
    print("M', dM:", M_prime, delta_m)
    output = open(ResultsName, 'w')
    for i in range(len(alpha_save)):
        output.write(str(alpha_save[i])+ " " + str(beta_save[i]) \
            + " " + str(Mp_save[i]) + " " + str(dM_save[i]) + " "\
            + str(om_save[i]) + " " + str(ol_save[i])+ " " + str(l_save[i])+ " " + str(w_save[i])+"\n")
    output.close()
    #print(alpha_save,beta_save,Mp_save,dM_save,om_save,ol_save)

if __name__ == "__main__":
    main()
