#!/usr/bin/python


# data wrapper for LC data from arXiv:1401.4064v2

import csv
import scipy as sc
import numpy as np

h = 0.673

H0 = 100 * h #km/s/Mpc


def inverseHubble(z,lm):
    # Calculates reduced hubble 
    return np.sqrt(lm*np.power((1+z),3)+(1-lm))

def dLumi(zhel, zcmb, lm):
    # takes heliocentric redshift + CMB frame redshift values
    # and then outputs a luminosity distance!
    # arxiv 1601.01451 for reference, eq 2,3
    distance = (1+zhel)*c*sc.integrate.quad(func, 0, zcmb, args=(lm,))/H0
    return distance

def import_data(inputfile):
    # blah blah
    datafile = open(inputfile, 'r')
    numLines = 0
    
    # returned lists
    mb_list = []
    dL = []
    c = []
    x1_list = []

    for line in datafile:
        if numLines == 0:
            continue
        dataline = list(csv.reader(f, delimiter=" "))
        #name zcmb zhel dz mb dmb x1 dx1 color dcolor 3rdvar d3rdvar tmax dtmax cov_m_s cov_m_c cov_s_c set ra dec biascor
        name = dataline[0]
        zcmb = float(dataline[1])
        zhel = float(dataline[2])
        dz = float(dataline[3])
        mb = float(dataline[4])
        dmb = float(dataline[5])
        x1 = float(dataline[6])
        dx1 = float(dataline[7])
        color = float(dataline[8])
        dcolor = float(dataline[9])

        dL = dLumi(zcmb,zhel)
        # append lists
        mb_list.append(mb)
        dL.append(dlumi)
        x1_list.append(x1)
        c.append(color)
    # returns list
    return mb_list, dL, c, x1_list
