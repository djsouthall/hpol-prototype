#!/usr/bin/env python3
'''
This will leverage some scripts from the BEACON work I've done (specifically)
ones that are used for reading in field fox data.  It will the plot the measured S11s.

Sloppy coding below for quick plotting.
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

sys.path.append(os.environ['BEACON_ANALYSIS_DIR'])
import tools.field_fox as ff

import matplotlib.pyplot as plt
from matplotlib import lines

def gainFromS21(filename_logmag, distance_m,header=17):
    '''
    Calculates the gain from the S21 following the source:
    https://www.pasternack.com/t-calculator-fspl.aspx

    This assumes the two antennas have the same gain. 
    '''
    try:
        freqs_Hz, logmag_vals = ff.readerFieldFox(filename_logmag,header=header)
        fspl = numpy.abs(logmag_vals)
        p1 = 20*numpy.log10(distance_m)
        p2 = 20*numpy.log10(freqs_Hz)
        p3 = 20*numpy.log10(4*numpy.pi/(299792458))
        gain = (p1 + p2 + p3 - fspl)/2.0
        #import pdb; pdb.set_trace()

        return freqs_Hz, gain,
    except Exception as e:
            print(e)

def plotGain(distance_m):
    '''
    '''
    plt.close('all')
    datapath ='/home/dsouthall/Projects/Greenland/hpol-prototype/data/Feb27-28_s21/'
    infiles = numpy.array(glob.glob(datapath + '*LOGMAG*.csv'))

    plt.rcParams['axes.prop_cycle'].by_key()['color']

    all_linestyles = numpy.array(list(lines.lineStyles.keys()))
    linestyles = numpy.array(['-', '--',':', '-.', ' ', ''], dtype='<U4')#all_linestyles[~numpy.isin(all_linestyles,['None'])]


    #PLOT PREPPING
    fontsize=16
    leg_fontsize=14
    alpha = 0.8
    #PLOT Gain
    gain_plot = plt.figure()
    gain_ax = plt.subplot(1,1,1)
    plt.ylabel('dBi',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    for infile in infiles:
        try:

            label = infile.split('/')[-1].replace('.csv','').replace('_PHASE','').replace('_LOGMAG','').replace('HOURGLASS','HOURGLASS ').replace('_',' ').title().replace('Pf',' pF').replace('Nh',' nH').replace('Ohm',' $\\Omega$').replace('Noshield','No Shield')

            freqs, gain = gainFromS21(infile,distance_m,header=17)
            
            plot_cut_ll = 100            
            plot_cut = freqs/1e6 > plot_cut_ll

            if 'No Shield' in label:
                linestyle = '--'
            else:
                linestyle = '-'


            gain_ax.plot(freqs[plot_cut]/1e6, gain[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)

        except Exception as e:
            print(e)


    gain_ax.legend(fontsize=leg_fontsize)
    gain_ax.set_xlim([plot_cut_ll,1000])

if __name__ == '__main__':
    plotGain(6.7818)