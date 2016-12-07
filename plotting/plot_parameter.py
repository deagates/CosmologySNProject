#!/usr/bin/python


# plotter for results code



import csv
#import scipy as sc
import numpy as np
import pylab as plt

final_om = []
final_weight = []
final_chi = []

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
            global final_om
            global final_weight
            global final_chi
            final_om.append(om)
            final_weight.append(w)
            final_chi.append(p)
            iterLines = 0
#            break

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
    main("results/Dec6_MCMC_Run_01.txt")
    main("results/Dec6_MCMC_Run_02.txt")
    main("results/Dec6_MCMC_Run_03.txt")
    main("results/Dec6_MCMC_Run_04.txt")
    main("results/Dec6_MCMC_Run_05.txt")
    main("results/Dec6_MCMC_Run_06.txt")
    main("results/Dec6_MCMC_Run_07.txt")
    main("results/Dec6_MCMC_Run_08.txt")
    main("results/Dec6_MCMC_Run_09.txt")
    main("results/Dec6_MCMC_Run_10.txt")
    print final_om, final_weight, final_chi
    weighted_omega_m_avg = sum([final_om[i]*final_weight[i] for i in range(len(final_om))])
    total_weight = sum([final_weight[i] for i in range(len(final_om))])
    print "WEIGHTED OMEGA_M AVG: ", weighted_omega_m_avg/total_weight
    plt.rcParams.update({'font.size': 18})
    n, bins, patches = plt.hist(final_om,10,normed=1,facecolor='pink')
    plt.xlabel(r"$\Omega_m$")
    plt.ylabel("A.U.")
    plt.savefig("plots/om_hist.pdf")
    plt.clf()

    plt.rcParams.update({'font.size': 18})
    plt.plot(final_om,final_chi,"b*")
    plt.xlabel(r"$\Omega_m$")
    plt.ylabel(r"$\chi^2$")
    #    plt.axis([0,n,0.,1.])
    plt.savefig("plots/om_chi.pdf")
#    main("results/Dec6_MCMC_Run_08.txt")

