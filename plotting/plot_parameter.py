#!/usr/bin/python


# plotter for results code



import csv
#import scipy as sc
import numpy as np
import pylab as plt
from scipy.stats import norm
import matplotlib.mlab as mlab

final_om = []
final_weight = []
final_chi = []
final_alpha = []
final_beta = []
final_dM = []
final_Mp = []
iters = 0
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

    om_save = 0
    ol_save = 0
    p_save = 0
    w_save = 0
    dM_save = 0
    Mp_save = 0
    


    for line in datafile:
        numLines = numLines + 1
        if numLines < 4:
            continue
        iterLines = iterLines + 1
        if line == "new iteration\n":
            iterLines = 0
            global iters
            iters = iters + 1
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
        if iterLines > 9999:
            global final_om
            global final_weight
            global final_chi
            global final_alpha
            global final_beta
            global final_dM
            global final_Mp
            if om_save == om and len(final_om) > 0:
                # if we had this value already, lets remove it
                final_om.pop()
                final_weight.pop()
                final_chi.pop()
                final_alpha.pop()
                final_beta.pop()
                final_dM.pop()
                final_Mp.pop()
            final_om.append(om)
            final_weight.append(w)
            final_chi.append(p)
            final_alpha.append(alpha)
            final_beta.append(beta)
            final_dM.append(dM)
            final_Mp.append(Mp)
        if iterLines == 20000:
            iterLines = 0
            #break
        om_save = om
        ol_save = ol
        p_save = p
        w_save = w
        dM_save = dM
        Mp_save = Mp
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
    #    print final_om, final_weight, final_chi
    print len(final_om), iters
    print "omega_m!"
    om_array = np.array(final_om)
    print np.percentile(om_array,50)
    print np.percentile(om_array,84)-np.percentile(om_array,50)
    print np.percentile(om_array,50)-np.percentile(om_array,16)

    print "alpha!"
    alpha_array = np.array(final_alpha)
    print np.percentile(alpha_array,50)
    print np.percentile(alpha_array,84)-np.percentile(alpha_array,50)
    print np.percentile(alpha_array,50)-np.percentile(alpha_array,16)

    print "beta!"
    beta_array = np.array(final_beta)
    print np.percentile(beta_array,50)
    print np.percentile(beta_array,84)-np.percentile(beta_array,50)
    print np.percentile(beta_array,50)-np.percentile(beta_array,16)

    print "dM!"
    dM_array = np.array(final_dM)
    print np.percentile(dM_array,50)
    print np.percentile(dM_array,84)-np.percentile(dM_array,50)
    print np.percentile(dM_array,50)-np.percentile(dM_array,16)

    print "Mp!"
    Mp_array = np.array(final_Mp)
    print np.percentile(Mp_array,50)
    print np.percentile(Mp_array,84)-np.percentile(Mp_array,50)
    print np.percentile(Mp_array,50)-np.percentile(Mp_array,16)

    weighted_omega_m_avg = sum([final_om[i] for i in range(len(final_om))])
    weighted_alpha = sum([final_alpha[i] for i in range(len(final_om))])
    weighted_beta = sum([final_beta[i] for i in range(len(final_om))])
    weighted_dM = sum([final_dM[i] for i in range(len(final_om))])
    weighted_Mp = sum([final_Mp[i] for i in range(len(final_om))])
    total_weight = len(final_Mp)
    # weighted_omega_m_avg = sum([final_om[i]*final_weight[i] for i in range(len(final_om))])
    # weighted_alpha = sum([final_alpha[i]*final_weight[i] for i in range(len(final_om))])
    # weighted_beta = sum([final_beta[i]*final_weight[i] for i in range(len(final_om))])
    # weighted_dM = sum([final_dM[i]*final_weight[i] for i in range(len(final_om))])
    # weighted_Mp = sum([final_Mp[i]*final_weight[i] for i in range(len(final_om))])
    # total_weight = sum([final_weight[i] for i in range(len(final_om))])
    print "WEIGHTED OMEGA_M AVG: ", weighted_omega_m_avg/total_weight
    print "WEIGHTED alpha AVG: ", weighted_alpha/total_weight
    print "WEIGHTED beta AVG: ", weighted_beta/total_weight
    print "WEIGHTED dM AVG: ", weighted_dM/total_weight
    print "WEIGHTED Mp AVG: ", weighted_Mp/total_weight
    plt.rcParams.update({'font.size': 18})
    n, bins, patches = plt.hist(final_om,100,normed=1,facecolor='pink')
    plt.xlabel(r"$\Omega_m$")
    plt.ylabel("A.U.")
    plt.savefig("plots/om_hist.pdf")
    plt.clf()

    plt.rcParams.update({'font.size': 18})
    n, bins, patches = plt.hist(final_alpha,100,normed=1,facecolor='pink')
    plt.xlabel(r"$\alpha$")
    plt.ylabel("A.U.")
    plt.savefig("plots/alpha_hist.pdf")
    plt.clf()

    plt.rcParams.update({'font.size': 18})
    n, bins, patches = plt.hist(final_beta,100,normed=1,facecolor='pink')
    plt.xlabel(r"$\beta$")
    plt.ylabel("A.U.")
    plt.savefig("plots/beta_hist.pdf")
    plt.clf()

    plt.rcParams.update({'font.size': 18})
    n, bins, patches = plt.hist(final_dM,100,normed=1,facecolor='pink')
    plt.xlabel(r"$\Delta_M$")
    plt.ylabel("A.U.")
    plt.savefig("plots/dM_hist.pdf")
    plt.clf()

    plt.rcParams.update({'font.size': 18})
    n, bins, patches = plt.hist(final_Mp,100,normed=1,facecolor='pink')
    plt.xlabel(r"$M_B'$")
    plt.ylabel("A.U.")
    plt.savefig("plots/Mp_hist.pdf")
    plt.clf()

    plt.rcParams.update({'font.size': 18})
    plt.plot(final_om,final_weight,"b*")
    plt.xlabel(r"$\Omega_m$")
#    plt.ylabel(r"$\chi^2$")
    plt.ylabel("Weight")
    #    plt.axis([0,n,0.,1.])
    plt.savefig("plots/om_weight.pdf")
#    main("results/Dec6_MCMC_Run_08.txt")

