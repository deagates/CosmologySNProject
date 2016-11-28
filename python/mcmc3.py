import sys, data_wrapper_dummy, random, math
from scipy.integrate import quad as integral
#dsfsfsfsf
final_params_omega_l = []
final_params_omega_m = []
final_params_omega_k = []
final_l = []
final_weights = []

c = 1
H0 = 1
dH=c/H0
M = 1

x1_data
m_B_data
C_data
z_data

# for now treating alpha and beta as constants that do not have to be fit for
alpha=0
beta=0

data_length=len(z_data)

Single_run_length=10
NUM_ITER = 10

def likelihood_calc(mu_model,mu_data):
    # evaluates likelihood given the model

    #using ChiSq
    likelihood=[]
    P = 0
    
    n = range(data_length-1)
    for i in n:
        p= (mu_data-mu_model)*(mu_data-mu_model)/mu_model
        P=P+p;
    return P


def omega_draw():
    om=random.random()
    ol=1-om
    ok=1
    #ok=floor(rand*3)-1
    return om, ol, ok

def H(x,omega_m,omega_l,omega_k):
    h = H0*math.sqrt(omega_m*(1+x)**3+omega_k*(1+x)**2+omega_l);
    return h

def mu_model_calc(z,_data,omega_m,omega_l,omega_k):
    # calculates mu given z
    mu_model=[];
    
    def ii(x):
        return H0/H(x,omega_m,omega_l,omega_k)
    
    n = range(data_length-1)
    for i in n:
        z=z_data[i]
        I, err=integral(ii(x),0,z)
        dL=(1+z)*(dH/sqrt(omega_k))*math.asin(sqrt(omgea_k)*I)
        f_i=25+ 5*log(dL)
        mu_model.append(f)
    return mu_model

def mu_data_calc(m_B,x1,C):
    # gives a functional form of mu given data + nuisance parameters
    mu_data=[]
    
    n = range(data_length-1)
    for i in n:
        m_B_i=m_B[i]
        x1_i_i=x1[i]
        C_i=C[i]
        f=m_B_i-M+alpha*x1_i-beta*C_i
        mu_data.append(f)
    print mu_data



def mu_error:
    # gives error associated
    

def likelihood_draw(z_data,m_B_data,x1_data,c_data)
    om,ol,ok=omega_draw()
    mu=mu_model_calc(z_data,om,ok,ol)
    mu_hat=mu_data_calc(m_B_data,x1_data,c_data)
    P=likelihood(mu_data,mu_model)
    return P,mu, mu_hat,om,ol,ok


def step_mcmc(Single_run_length,z_data,m_B_data,x1_data,c_data):
    # evaluates one step of mcmc algorithm
    # put weight here for now
    #initialize 
    P_old,mu, mu_hat,om,ol,ok = likelihood_draw(z_data,m_B_data,x1_data,c_data);
    w=0
    n = range(Single_run_length)
    for i in n:
        P_new,mu, mu_hat,om,ol,ok = likelihood_draw(z_data,m_B_data,x1_data,c_data)
        r=P_new/P_old;
        if r>1:
            w=w+1;
        else:
            w=0;
            P_old=P_new;
    return om, ol, ok, w, P_new

    
def covariance:
    # think about this later

def main:
    # call functions
    n = range(NUM_ITER)
    for i in n:
        param_omega_m, param_omega_l, param_omega_k, weight, likelihood = step_mcmc(Single_run_length,z_data,m_B_data,x1_data,c_data)
        final_params_omega_m.append(param_omega_m)
        final_params_omega_l.append(param_omega_l)
        final_params_omega_k.append(param_omega_k)
        final_l.append(likelihood)
        final_weights.append(weight)
    print final_params_omega_m,final_params_omega_l,final_params_omega_k, final_l, final_weights
