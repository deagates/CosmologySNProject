#!/usr/bin/python


# plotter for results code



import csv
#import scipy as sc
import numpy as np
import pylab as plt


def import_data(inputfile):
    # expects: 
    # alpha, beta, delta_M, M_prime, omega_m, omega_l, likelihood
    datafile = open(inputfile, 'r')
    numLines = 0
    iterLines = 0
        # returned lists
    alpha_list = []
    beta_list = []
    dM_list = []
    Mp_list = []
    om_list = []
    ol_list = []
    p_list = []
    w_list = []

    final_om = []
    final_weight = []

    for line in datafile:
        numLines = numLines + 1
        if numLines < 4:
            continue
        iterLines = iterLines + 1
        if line == "new iteration\n":
            iterLines = 0
            continue
        dataline = line.split()
        alpha = float(dataline[0])
        beta = float(dataline[1])
        Mp = float(dataline[2])
        dM = float(dataline[3])
        om = float(dataline[4])
        ol = float(dataline[5])
        p = float(dataline[6])
        w = float(dataline[7])
        alpha_list.append(alpha)
        beta_list.append(beta)
        dM_list.append(dM)
        Mp_list.append(Mp)
        om_list.append(om)
        ol_list.append(ol)
        p_list.append(p)
        w_list.append(w)
        if iterLines == 20000:
            final_om.append(om)
            final_weight.append(w)
            iterLines = 0
#            break

    print final_weight, final_om
    
    return alpha_list,beta_list,dM_list,Mp_list,om_list,ol_list,p_list

def main(inputfile):
    a,b,dM,Mp,om,ol,p = import_data(inputfile)
    # n = len(a)
    # plt.rcParams.update({'font.size': 18})
    # plt.plot(range(n),a)
    # plt.xlabel("Steps")
    # plt.ylabel(r"$\alpha$")
    # plt.axis([0,n,-1.,4.])
    # plt.savefig("plots/alpha.pdf")    
    # plt.clf()

    # plt.plot(range(n),b)
    # plt.xlabel("Steps")
    # plt.ylabel(r"$\beta$")
    # plt.axis([0,n,1.,6.])
    # plt.savefig("plots/beta.pdf")
    # plt.clf()

    # plt.plot(range(n),dM)
    # plt.xlabel("Steps")
    # plt.ylabel(r"$\Delta$M")
    # plt.axis([0,n,-2.,1.])
    # plt.savefig("plots/dM.pdf")    
    # plt.clf()

    # plt.plot(range(n),Mp)
    # plt.xlabel("Steps")
    # plt.ylabel("M'")
    # plt.axis([0,n,-21.,-17.])
    # plt.savefig("plots/Mp.pdf")
    # plt.clf()

    # plt.plot(range(n),om)
    # plt.xlabel("Steps")
    # plt.ylabel(r"$\Omega_m$")
    # plt.axis([0,n,0.,1.])
    # plt.savefig("plots/om.pdf")
    # plt.clf()

    # plt.plot(range(n),ol)
    # plt.xlabel("Steps")
    # plt.ylabel(r"$\Omega_l$")
    # plt.axis([0,n,0.,1.])
    # plt.savefig("plots/ol.pdf")
if __name__ == "__main__":
    main("results/Dec6_MCMC_Run_06.txt")
#    main("results/Dec6_MCMC_Run_08.txt")

