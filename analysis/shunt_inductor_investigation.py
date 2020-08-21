#!/usr/bin/env python3
'''
This is meant to be the code that searches for the best shunt inductor value to try for the hpol prototype.
It will may also go down the road of investigating discrepancy between measured values in the lab and those
predicted by the Smith chart calculations. 
'''

import numpy
import scipy.spatial
import scipy.signal
from scipy.optimize import curve_fit
import os
import sys
import csv
import inspect
import glob

from beacon.tools import field_fox as ff
from hpol_prototype.tools.matching_networks import HpolTriWingFeed as Feed

import matplotlib.pyplot as plt
from matplotlib import lines



if __name__ == '__main__':
    plt.close('all')
    datapath ='/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_copper_tab_testing_aug2020/s11/'
    infiles = numpy.array(glob.glob(datapath + '*.csv'))
    roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])


    #Prepare the starting S11 that is to be shifted around.
    for root in roots:
        if 'NOSHUNT' in root:
            starter_root = root
    roots = roots[roots != starter_root]
    starter_logmag_infile = starter_root.replace('_FILLER','_LOGMAG')
    starter_phase_infile = starter_root.replace('_FILLER','_PHASE')

    #Load intial data and make feed
    starter_freqs, starter_logmag_vals = ff.readerFieldFox(starter_logmag_infile,header=17)
    starter_freqs, starter_phase_vals = ff.readerFieldFox(starter_phase_infile,header=17)
    feed = Feed(starter_freqs, starter_logmag_vals, starter_phase_vals, z_0=50)

    #inductor_values_nH = numpy.linspace(15,250,300)
    inductor_values_nH = numpy.linspace(5.0,200.0,500)#numpy.array([47.0,56,100,110])
    test_statistic = numpy.zeros_like(inductor_values_nH)
    db_cut = -5 #dB

    for index, L in enumerate(inductor_values_nH):
        feed.addShuntInductor(L*1.0e-9)
        test_statistic[index] = feed.getPercentBelow(db=db_cut)
        feed.reset()


    plt.figure()
    plt.plot(inductor_values_nH, test_statistic,label='Predicted Curve from Varying No Shunt Feed')
    plt.xlabel('Individual Added Shunt Inductor Values (nH)')
    plt.ylabel('Fraction of S11 Points below %0.2f dB'%db_cut)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-',alpha=0.75)
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.25)
    feed.reset()
    plt.axhline(feed.getPercentBelow(db=db_cut),c='r',linestyle='--',label='Measured No Shunt Fraction')

    measured_test_statistic = numpy.zeros(len(roots))
    measured_inductor_values_nH = numpy.zeros(len(roots))
    for index, root in enumerate(roots):
        logmag_infile = root.replace('_FILLER','_LOGMAG')
        phase_infile = root.replace('_FILLER','_PHASE')
        L = root.split('/')[-1].split('_')[2]
        if L.lower() == 'noshunt':
            L = 0
        else:
            L = float(L.lower().replace('nh',''))
        measured_inductor_values_nH[index] = L
        freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
        freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)
        measured_feed = Feed(freqs, logmag_vals, phase_vals, z_0=50)
        measured_test_statistic[index] = measured_feed.getPercentBelow(db=db_cut)
        #measured_feed.plotCurrentLogMagS11(label=root)



    plt.scatter(measured_inductor_values_nH,measured_test_statistic,color='r',edgecolors='k',label='Measured Feed Values')
    plt.legend(loc='upper right')

    # print('Testing low L values')
        
    # test_L = 5e-9

    # test_feed = Feed(starter_freqs, starter_logmag_vals, starter_phase_vals, z_0=50)

    # test_feed.reset()
    # ax1 = test_feed.plotCurrentLogMagS11(label='Original')
    # test_feed.addShuntInductor(test_L)
    # test_feed.plotCurrentLogMagS11(ax=ax1,label='Adjusted Test')

    # test_feed.reset()
    # ax2 = test_feed.plotSmithChart(label='Original')
    # test_feed.addShuntInductor(test_L)
    # test_feed.plotSmithChart(ax=ax2,label='Adjusted Test')


    #testing per wing adding v.s. equivalent math

    test_feed = Feed(starter_freqs, starter_logmag_vals, starter_phase_vals, z_0=50)
    ax = test_feed.plotCurrentLogMagS11(label='original')
    test_feed.addRLC('pl',45e-9)
    test_feed.addRLC('pl',45e-9)
    test_feed.addRLC('pl',45e-9)
    ax = test_feed.plotCurrentLogMagS11(ax=ax,label='added one at a time')
    test_feed.reset()
    test_feed.addRLC('pl',(45e-9)/3)
    ax = test_feed.plotCurrentLogMagS11(ax=ax,label='added equivalent L')

