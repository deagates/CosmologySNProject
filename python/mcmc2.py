import sys, data_wrapper_dummy
#dsfsfsfsf
final_params_omega_l = []
final_params_omega_m = []
final_l = []
final_weights = []

c = 1
H0 = 1
dH=c/H0
M = 1

x1
m_B
C
z_data

data_length=len(z_data)


NUM_ITER = 10

def likelihood_calc(mu_model,mu_data):
    # evaluates likelihood given the model

    #using ChiSq
    likelihood=[];
for i in data_length
#chi2
P=0;

for i=1:data_length
    p= (mu_data-mu_model)^2/mu_model;
    P=P+p;
end
end



def omega_draw:

om=random()
ol=1-om
ok=1
#ok=floor(rand*3)-1

print om, ol, ok

def H(omega_m,omega_l,omega_k):
    h = H0*sqrt(omega_m*(1+x)^3+omega_k*(1+x)^2+omega_l);
    print h

def mu_model_calc(z,_data,omega_m,omega_l,omega_k):
    # calculates mu given z
mu_model=[];
for i in data_length
    z=z_data[i]
    I=integral(H0/H(z,omega_m,omega_l,omega_k),0,z)
    dL=(1+z)*(dH/sqrt(omega_k))*asin(sqrt(omgea_k)*I)
    f_i=25+ 5*log(dL)
    mu_model.append(f)
end
print mu_model

def mu_data_calc(m_B,x1,C):
    # gives a functional form of mu given data + nuisance parameters
mu_data=[]

n = range(data_length)
for i in data_length
    m_B_i=m_B[i]
    x1_i_i=x1[i]
    C_i=C[i]
    f=m_B_i-M+alpha*x1_i-beta*C_i
    mu_data.append(f)
    print mu_data
end


def mu_error:
    # gives error associated
    

def likelihood_draw(z_data,m_B_data,x1_data,c_data)
om,ol,ok=omega_draw();
mu=mu_model_calc(z_data,om,ok,ol);
mu_hat=mu_data_calc(m_B_data,x1_data,c_data)
P=likelihood(mu_data,mu_model);

print P,mu, mu_hat,om,ol,ok


def step_mcmc(Single_run_length,z_data,m_B_data,x1_data,c_data):
    # evaluates one step of mcmc algorithm
    # put weight here for now
%initialize 
[P_old,mu, mu_hat,om,ol,ok] = likelihood_draw(z_data,m_B_data,x1_data,c_data);
w=0;

for ii=1:Single_run_length
    
    [P_new,mu, mu_hat,om,ol,ok] = likelihood_draw(z_data,m_B_data,x1_data,c_data)
    
    r=P_new/P_old;
    
    if r<1
        w=w+1;
    else
        w=0;
        P_old=P+new;
    end
end

olfinal=ol;
omfinal=om;



    # return p, l, w
    
def covariance:
    # think about this later

def main:
    # call functions
    n = range(NUM_ITER)
    for i in n:
        param, likelihood, weight = step_mcmc()
        final_params.append(param)
        final_l.append(likelihood)
        final_weights.append(weight)
    print final_params, final_l, final_weights
