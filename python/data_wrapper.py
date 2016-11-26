#!/usr/bin/python


# data wrapper for LC data from arXiv:1401.4064v2

import csv
import scipy as sc
import numpy as np

h = 0.673

H0 = 100 * h #km/s/Mpc


def inverseHubble(z,lm):
    # Calculates reduced hubble assuming lambda CDM
    # lm is omega_matter
    return np.sqrt(lm*np.power((1+z),3)+(1-lm))

def abs_mag(hostmass):
    # returns true or false given hostmass
    # arXiv:1401.4064v2 proposed step function for handling differences
    # for various host galaxies
    # useless function but this is just to remind me
    # to put this in the calculation for mu
    if (hostmass < 10.):
        return true
    else:
        return false

def dLumi(zhel, zcmb, lm, func):
    # takes heliocentric redshift + CMB frame redshift values
    # and then outputs a luminosity distance!
    # arxiv 1601.01451 for reference, eq 2,3
    distance = (1+zhel)*c*sc.integrate.quad(func, 0, zcmb, args=(lm,))/H0
    return distance

def import_data(inputfile):
    datafile = open(inputfile, 'r')
    numLines = 0
    
    # returned lists
    mb_list = []
    zcmb_list = []
    zhel_list = []
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
        hostmass = float(dataline[10])
        dhostmass = float(dataline[11])

#        dL = dLumi(zcmb,zhel)

# calculating dL requires a guess for omega_m, so I can't do it here, should be called in mcmc,
# but you can try using dLumi, inverseHubble

        # append lists
        mb_list.append(mb)
        zcmb_list.append(zcmb)
        zhel_list.append(zhel)
        x1_list.append(x1)
        c.append(color)
    # returns list
    return mb_list, zcmb, zhel, c, x1_list

