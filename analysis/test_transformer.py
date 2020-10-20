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
from hpol_prototype.tools.matching_networks import SmithMatcher
from hpol_prototype.tools.pySmith import get_smith

import matplotlib.pyplot as plt
from matplotlib import lines


def makeGeneralPlots3():
    '''
    This will make s11 expected v.s. measured plots in both logmag and smith chart.   As well as the general
    curve for expected points below threshold.  
    '''
    


if __name__ == '__main__':
    readout_transformer_factor = 1.5
    plot_transformed = True
    try:
        if False:
            datapath ='/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_al_tube_testing_sep2020/s11/'
            starter_root_match = 'hpol_ernieemu_100nh_2p7pf'
        elif False:
            datapath ='/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_al_tube_testing_sep2020/s11/'
            starter_root_match = 'hpol_ernieemu_22nh_2p7pF'
        elif False:
            datapath ='/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_copper_tab_testing_aug2020/s11/'
            starter_root_match = 'hpol_ape_110nh_2p7pf'
        elif False:
            datapath ='/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_copper_tab_testing_aug2020/s11/'
            starter_root_match = 'hpol_bee_56nh_2p7pf'
        elif False:
            datapath ='/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_copper_tab_testing_aug2020/s11/'
            starter_root_match = 'hpol_bee_47nh_2p7pf'
        elif True:
            datapath ='/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_copper_tab_testing_aug2020/s11/'
            starter_root_match = 'hpol_chuckape_27nh_2p7pf'
        elif False:
            datapath ='/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_copper_tab_testing_aug2020/s11/'
            starter_root_match = 'hpol_davebee_27nh_2p7pf'
        elif False:
            datapath ='/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_copper_tab_testing_aug2020/s11/'
            starter_root_match = 'hpol_bee_22nh_2p7pf'
        else:
            datapath ='/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_copper_tab_testing_aug2020/s11/'
            starter_root_match = 'hpol_bee_NOSHUNT_2p7pF'

        #starter_root_match = datapath + starter_root_match + '_FILLER.csv'

        #inductor_values_nH = numpy.linspace(15,250,300)
        fontsize=16

        #Prep calculations
        starter_shunt_inductor_value = starter_root_match.split('/')[-1].split('_')[2] #nH
        if starter_shunt_inductor_value.lower() == 'noshunt':
            starter_shunt_inductor_value = 0
        else:
            starter_shunt_inductor_value = float(starter_shunt_inductor_value.lower().replace('nh',''))*1e-9 #H
        starter_cap_value = float(starter_root_match.split('/')[-1].split('_')[3].lower().replace('pf','').replace('p','.'))*1e-12 #F

        #datapath ='/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_al_tube_testing_sep2020/s11/'
        infiles = numpy.array(glob.glob(datapath + '*.csv'))
        roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])
        
        #Prepare the starting S11 that is to be shifted around.
        for root in roots:
            if starter_root_match in root:
                starter_root = root
        roots = roots[roots != starter_root]
        starter_logmag_infile = starter_root.replace('_FILLER','_LOGMAG')
        starter_phase_infile = starter_root.replace('_FILLER','_PHASE')

        #Load intial data and make feed
        starter_freqs, starter_logmag_vals = ff.readerFieldFox(starter_logmag_infile,header=17)
        starter_freqs, starter_phase_vals = ff.readerFieldFox(starter_phase_infile,header=17)
        feed = SmithMatcher(starter_freqs, starter_logmag_vals, starter_phase_vals, initial_capacitor_value=starter_cap_value, initial_shunt_inductor_value=starter_shunt_inductor_value, z_0=50, readout_transformer_factor=readout_transformer_factor)

        # feed.plotCurrentLogMagS11(label='Untransformed',plot_transformed=False)
        # feed.plotCurrentLogMagS11(label='transformed',plot_transformed=True)

        compare_fig = plt.figure()
        compare_ax1 = plt.subplot(2,1,1)

        compare_ax1 = feed.plotCurrentLogMagS11(ax=compare_ax1,label='Untransformed',plot_transformed=False) #Labelling assumes starter_ curves are for no shunt curve.
        compare_ax1 = feed.plotCurrentLogMagS11(ax=compare_ax1,label='transformed',plot_transformed=True)
        plt.ylabel('dB',fontsize=14)
        plt.xlabel('MHz',fontsize=14)
        plt.minorticks_on()
        plt.grid(b=True, which='major', color='k', linestyle='-')
        plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

        feed.reset()
        #compare_ax2 = plt.subplot(2,1,2)
        compare_ax2 = get_smith(compare_fig, rect = 212, plot_impedance = False, plot_ticks = True, plot_admittance = True, plot_labels = False)
        plt.ylabel('Im($\\Gamma$)',fontsize=14)
        plt.xlabel('Re($\\Gamma$)',fontsize=14)
        compare_ax2 = feed.plotSmithChart(ax=compare_ax2,plot_transformed=False) #Labelling assumes starter_ curves are for no shunt curve.
        compare_ax2 = feed.plotSmithChart(ax=compare_ax2,plot_transformed=True)


    except Exception as e:
        print('\nError in %s'%inspect.stack()[0][3])
        print(e)
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)