import sys, data_wrapper
#dsfsfsfsf
final_params = []
final_l = []
final_weights = []
NUM_ITER = 10
def likelihood:
    # evaluates likelihood given the model

def mu_model:
    # calculates mu given z

def mu_data:
    # gives a functional form of mu given data + nuisance parameters

def mu_error:
    # gives error associated
    
def step_mcmc:
    # evaluates one step of mcmc algorithm
    # put weight here for now

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
