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
    
        # returned lists
    alpha_list = []
    beta_list = []
    dM_list = []
    Mp_list = []
    om_list = []
    ol_list = []
    p_list = []

    for line in datafile:
        numLines = numLines + 1
        if numLines < 5001:
            continue
        dataline = line.split()
        alpha = float(dataline[0])
        beta = float(dataline[1])
        dM = float(dataline[2])
        Mp = float(dataline[3])
        om = float(dataline[4])
        ol = float(dataline[5])
        p = float(dataline[6])

        alpha_list.append(alpha)
        beta_list.append(beta)
        dM_list.append(dM)
        Mp_list.append(Mp)
        om_list.append(om)
        ol_list.append(ol)
        p_list.append(p)

        if numLines == 10000:
            break

        # returns list
    return alpha_list,beta_list,dM_list,Mp_list,om_list,ol_list,p_list

def main(inputfile):
    a,b,dM,Mp,om,ol,p = import_data(inputfile)
    n = len(a)
    plt.plot(range(n),a)
    plt.xlabel("Steps")
    plt.ylabel("Alpha")
    plt.axis([0,n,0.,4.])
    plt.savefig("alpha.pdf")    
    plt.clf()

    plt.plot(range(n),b)
    plt.xlabel("Steps")
    plt.ylabel("Beta")
    plt.axis([0,n,1.,5.])
    plt.savefig("beta.pdf")
    plt.clf()

    plt.plot(range(n),dM)
    plt.xlabel("Steps")
    plt.ylabel(r"$\Delta$M")
    plt.axis([0,n,-1.,0.])
    plt.savefig("dM.pdf")    
    plt.clf()

    plt.plot(range(n),Mp)
    plt.xlabel("Steps")
    plt.ylabel("M'")
    plt.axis([0,n,-20.,-18.])
    plt.savefig("Mp.pdf")
    plt.clf()

    plt.plot(range(n),om)
    plt.xlabel("Steps")
    plt.ylabel(r"$\Omega_m$")
    plt.axis([0,n,0.,1.])
    plt.savefig("om.pdf")
    plt.clf()

    plt.plot(range(n),ol)
    plt.xlabel("Steps")
    plt.ylabel(r"$\Omega_l$")
    plt.axis([0,n,0.,1.])
    plt.savefig("ol.pdf")
if __name__ == "__main__":
    main("results3.txt")

