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


def test_abcd():
    '''
    '''
    # plt.close('all')
    datapath ='/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data_ABCD/'
    infiles = numpy.array(glob.glob(datapath + '*CAT*.csv'))
    infiles = numpy.append(infiles,numpy.array(glob.glob(datapath + '*DOG*.csv')))

    #caps = numpy.unique([infile.split('/')[-1].split('_')[1] for infile in infiles])
    #shapes = numpy.unique([infile.split('/')[-1].split('_')[2] for infile in infiles])

    plt.rcParams['axes.prop_cycle'].by_key()['color']

    all_linestyles = numpy.array(list(lines.lineStyles.keys()))
    linestyles = numpy.array(['-', '--',':', '-.', ' ', ''], dtype='<U4')#all_linestyles[~numpy.isin(all_linestyles,['None'])]
    roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])


    #PLOT PREPPING
    fontsize=16
    leg_fontsize=14
    #PLOT 1, LOGMAG
    logmag_plot = plt.figure()
    logmag_ax = plt.subplot(1,1,1)
    plt.ylabel('dB',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    #PLOT 2, Impedance
    impedance_plot = plt.figure()
    impedance_ax1 = plt.subplot(3,1,1)

    plt.ylabel('Re(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax2 = plt.subplot(3,1,2,sharex=impedance_ax1)

    plt.ylabel('Im(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax3 = plt.subplot(3,1,3,sharex=impedance_ax1)

    plt.ylabel('|Z|',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    ignore_roots = ['/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_3.3PF-NO-SHUNT_STRAIGHT_FILLER.csv',\
                    '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF_STRAIGHT_16-75IN.csv',\
                    '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF_STRAIGHT_15-30IN.csv',\
                    '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF_STRAIGHT_135LOGMAG.csv']

    ignore_roots = []

    for root in roots:
        try:
            if root in ignore_roots:
                continue
            print(root)
            if root == '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF-45NH-SHUNT_STRAIGHT_FILLER.csv':
                alpha=1.0
            else:
                alpha=0.8
            logmag_infile = root.replace('_FILLER','_LOGMAG')
            phase_infile = root.replace('_FILLER','_PHASE')

            label = logmag_infile.split('/')[-1].replace('.csv','').replace('_PHASE','').replace('_LOGMAG','').replace('HOURGLASS','HOURGLASS ').replace('_',' ').title().replace('Pf',' pF').replace('Nh',' nH').replace('Ohm',' $\\Omega$')
            #cap = logmag_infile.split('/')[-1].split('_')[1]
            #shape = logmag_infile.split('/')[-1].split('_')[2]

            #color = plt.rcParams['axes.prop_cycle'].by_key()['color'][numpy.where(numpy.isin(caps,[cap]))[0][0]]
            #linestyle = linestyles[numpy.where(numpy.isin(shapes,[shape]))[0][0]]

            freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
            freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)
            
            plot_cut_ll = 100            
            plot_cut = freqs/1e6 > plot_cut_ll

            if 'Cat' in label:
                linestyle = '--'
            elif 'Dog' in label:
                linestyle = '-'
            else:
                linestyle = None


            logmag_ax.plot(freqs[plot_cut]/1e6, logmag_vals[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)

            #data = ff.plotLogMag(logmag_infile,freq_range=(min(freqs),max(freqs)),fig=logmag_plot,label=label,color=color,linestyle=linestyle,alpha=0.8)


            lin_vals = ff.logMagToLin(logmag_vals)

            Re, Im = ff.magPhaseToReIm(lin_vals,phase_vals)

            complex_S11 = Re + 1j*Im

            z = ff.linToComplexZ( complex_S11, z_0 = 50.0 )

            impedance_ax1.plot(freqs[plot_cut]/1e6, numpy.real(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax2.plot(freqs[plot_cut]/1e6, numpy.imag(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax3.plot(freqs[plot_cut]/1e6, numpy.abs(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
        except Exception as e:
            print(e)





    logmag_ax.legend(fontsize=leg_fontsize)
    logmag_ax.set_xlim([plot_cut_ll,1000])


    impedance_ax1.legend(fontsize=leg_fontsize,loc='upper right')


    impedance_ax2.legend(fontsize=leg_fontsize,loc='upper right')

    impedance_ax3.legend(fontsize=leg_fontsize,loc='upper right')


def test_switch_to_consistent_wire_lengths():
    '''
    '''
    # plt.close('all')
    datapath ='/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data_ABCD/'
    infiles = numpy.array(glob.glob(datapath + '*2PF*ALICE-CAT*.csv'))
    infiles = numpy.append(infiles,numpy.array(glob.glob(datapath + '*2PF*BOB-DOG*.csv')))

    #caps = numpy.unique([infile.split('/')[-1].split('_')[1] for infile in infiles])
    #shapes = numpy.unique([infile.split('/')[-1].split('_')[2] for infile in infiles])

    plt.rcParams['axes.prop_cycle'].by_key()['color']

    all_linestyles = numpy.array(list(lines.lineStyles.keys()))
    linestyles = numpy.array(['-', '--',':', '-.', ' ', ''], dtype='<U4')#all_linestyles[~numpy.isin(all_linestyles,['None'])]
    roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])


    #PLOT PREPPING
    fontsize=16
    leg_fontsize=14
    #PLOT 1, LOGMAG
    logmag_plot = plt.figure()
    logmag_ax = plt.subplot(1,1,1)
    plt.ylabel('dB',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    #PLOT 2, Impedance
    impedance_plot = plt.figure()
    impedance_ax1 = plt.subplot(3,1,1)

    plt.ylabel('Re(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax2 = plt.subplot(3,1,2,sharex=impedance_ax1)

    plt.ylabel('Im(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax3 = plt.subplot(3,1,3,sharex=impedance_ax1)

    plt.ylabel('|Z|',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    ignore_roots = ['/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_3.3PF-NO-SHUNT_STRAIGHT_FILLER.csv',\
                    '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF_STRAIGHT_16-75IN.csv',\
                    '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF_STRAIGHT_15-30IN.csv',\
                    '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF_STRAIGHT_135LOGMAG.csv']

    ignore_roots = []

    for root in roots:
        try:
            if root in ignore_roots:
                continue
            print(root)
            if root == '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF-45NH-SHUNT_STRAIGHT_FILLER.csv':
                alpha=1.0
            else:
                alpha=0.8
            logmag_infile = root.replace('_FILLER','_LOGMAG')
            phase_infile = root.replace('_FILLER','_PHASE')

            label = logmag_infile.split('/')[-1].replace('.csv','').replace('_PHASE','').replace('_LOGMAG','').replace('HOURGLASS','HOURGLASS ').replace('_',' ').title().replace('Pf',' pF').replace('Nh',' nH').replace('Ohm',' $\\Omega$')
            #cap = logmag_infile.split('/')[-1].split('_')[1]
            #shape = logmag_infile.split('/')[-1].split('_')[2]

            #color = plt.rcParams['axes.prop_cycle'].by_key()['color'][numpy.where(numpy.isin(caps,[cap]))[0][0]]
            #linestyle = linestyles[numpy.where(numpy.isin(shapes,[shape]))[0][0]]

            freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
            freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)
            
            plot_cut_ll = 100            
            plot_cut = freqs/1e6 > plot_cut_ll

            if 'Cat' in label:
                linestyle = '--'
            elif 'Dog' in label:
                linestyle = '-'
            else:
                linestyle = None


            logmag_ax.plot(freqs[plot_cut]/1e6, logmag_vals[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)

            #data = ff.plotLogMag(logmag_infile,freq_range=(min(freqs),max(freqs)),fig=logmag_plot,label=label,color=color,linestyle=linestyle,alpha=0.8)


            lin_vals = ff.logMagToLin(logmag_vals)

            Re, Im = ff.magPhaseToReIm(lin_vals,phase_vals)

            complex_S11 = Re + 1j*Im

            z = ff.linToComplexZ( complex_S11, z_0 = 50.0 )

            impedance_ax1.plot(freqs[plot_cut]/1e6, numpy.real(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax2.plot(freqs[plot_cut]/1e6, numpy.imag(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax3.plot(freqs[plot_cut]/1e6, numpy.abs(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
        except Exception as e:
            print(e)





    logmag_ax.legend(fontsize=leg_fontsize)
    logmag_ax.set_xlim([plot_cut_ll,1000])


    impedance_ax1.legend(fontsize=leg_fontsize,loc='upper right')


    impedance_ax2.legend(fontsize=leg_fontsize,loc='upper right')

    impedance_ax3.legend(fontsize=leg_fontsize,loc='upper right')


def test_47():
    '''
    '''
    # plt.close('all')
    datapath ='/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data_ABCD/'
    print(datapath)
    infiles = numpy.array(glob.glob(datapath + '*STRAIGHT*.csv'))
    infiles = numpy.append(infiles,numpy.array(glob.glob(datapath + '*DOG2*.csv')))
    infiles = numpy.append(infiles,numpy.array(glob.glob(datapath + '*CAT2*.csv')))
    print(infiles)

    #caps = numpy.unique([infile.split('/')[-1].split('_')[1] for infile in infiles])
    #shapes = numpy.unique([infile.split('/')[-1].split('_')[2] for infile in infiles])

    plt.rcParams['axes.prop_cycle'].by_key()['color']

    all_linestyles = numpy.array(list(lines.lineStyles.keys()))
    linestyles = numpy.array(['-', '--',':', '-.', ' ', ''], dtype='<U4')#all_linestyles[~numpy.isin(all_linestyles,['None'])]
    roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])


    #PLOT PREPPING
    fontsize=16
    leg_fontsize=14
    #PLOT 1, LOGMAG
    logmag_plot = plt.figure()
    logmag_ax = plt.subplot(1,1,1)
    plt.ylabel('dB',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    #PLOT 2, Impedance
    impedance_plot = plt.figure()
    impedance_ax1 = plt.subplot(3,1,1)

    plt.ylabel('Re(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax2 = plt.subplot(3,1,2,sharex=impedance_ax1)

    plt.ylabel('Im(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax3 = plt.subplot(3,1,3,sharex=impedance_ax1)

    plt.ylabel('|Z|',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    ignore_roots = ['/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_3.3PF-NO-SHUNT_STRAIGHT_FILLER.csv',\
                    '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF_STRAIGHT_16-75IN.csv',\
                    '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF_STRAIGHT_15-30IN.csv',\
                    '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF_STRAIGHT_135LOGMAG.csv']

    ignore_roots = []

    for root in roots:
        try:
            if root in ignore_roots:
                continue
            print(root)
            if root == '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF-45NH-SHUNT_STRAIGHT_FILLER.csv':
                alpha=1.0
            else:
                alpha=0.8
            logmag_infile = root.replace('_FILLER','_LOGMAG')
            phase_infile = root.replace('_FILLER','_PHASE')

            label = logmag_infile.split('/')[-1].replace('.csv','').replace('_PHASE','').replace('_LOGMAG','').replace('HOURGLASS','HOURGLASS ').replace('_',' ').title().replace('Pf',' pF').replace('Nh',' nH').replace('Ohm',' $\\Omega$')
            #cap = logmag_infile.split('/')[-1].split('_')[1]
            #shape = logmag_infile.split('/')[-1].split('_')[2]

            #color = plt.rcParams['axes.prop_cycle'].by_key()['color'][numpy.where(numpy.isin(caps,[cap]))[0][0]]
            #linestyle = linestyles[numpy.where(numpy.isin(shapes,[shape]))[0][0]]

            freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
            freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)
            
            plot_cut_ll = 100            
            plot_cut = freqs/1e6 > plot_cut_ll

            if 'Cat' in label:
                linestyle = '--'
            elif 'Dog' in label:
                linestyle = '-'
            else:
                linestyle = None


            logmag_ax.plot(freqs[plot_cut]/1e6, logmag_vals[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)

            #data = ff.plotLogMag(logmag_infile,freq_range=(min(freqs),max(freqs)),fig=logmag_plot,label=label,color=color,linestyle=linestyle,alpha=0.8)


            lin_vals = ff.logMagToLin(logmag_vals)

            Re, Im = ff.magPhaseToReIm(lin_vals,phase_vals)

            complex_S11 = Re + 1j*Im

            z = ff.linToComplexZ( complex_S11, z_0 = 50.0 )

            impedance_ax1.plot(freqs[plot_cut]/1e6, numpy.real(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax2.plot(freqs[plot_cut]/1e6, numpy.imag(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax3.plot(freqs[plot_cut]/1e6, numpy.abs(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
        except Exception as e:
            print(e)
    logmag_ax.legend(fontsize=leg_fontsize)
    logmag_ax.set_xlim([plot_cut_ll,1000])


    impedance_ax1.legend(fontsize=leg_fontsize,loc='upper right')


    impedance_ax2.legend(fontsize=leg_fontsize,loc='upper right')

    impedance_ax3.legend(fontsize=leg_fontsize,loc='upper right')

def test_near45(mode=1):
    '''
    '''
    datapath ='/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data_near45/'
    print(datapath)
    plt.rcParams['axes.prop_cycle'].by_key()['color']

    all_linestyles = numpy.array(list(lines.lineStyles.keys()))
    linestyles = numpy.array(['-', '--',':', '-.', ' ', ''], dtype='<U4')#all_linestyles[~numpy.isin(all_linestyles,['None'])]
    if mode == 0:
        # plt.close('all')
        infiles = numpy.array(glob.glob(datapath + '*47*CAT*.csv'))
        infiles = numpy.append(infiles,numpy.array(glob.glob(datapath + '*33*CAT*.csv')))
        infiles = numpy.sort(infiles)
        roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])
        infiles = numpy.append(infiles,numpy.array(glob.glob(datapath + '*45*.csv')))
        _roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in numpy.array(glob.glob(datapath + '*45*.csv'))])
        print(roots)
        print(_roots)
        roots = numpy.append(roots,_roots)
        print(roots)
    elif mode == 1:
        infiles = numpy.array(glob.glob(datapath + '*47*CAT*.csv'))
        infiles = numpy.append(infiles,numpy.array(glob.glob(datapath + '*33*CAT*.csv')))
        #infiles = numpy.append(infiles,numpy.array(glob.glob(datapath + '*DOG2*.csv')))
        #infiles = numpy.append(infiles,numpy.array(glob.glob(datapath + '*CAT2*.csv')))
        infiles = numpy.sort(infiles)
        roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])
    print(infiles)

    #caps = numpy.unique([infile.split('/')[-1].split('_')[1] for infile in infiles])
    #shapes = numpy.unique([infile.split('/')[-1].split('_')[2] for infile in infiles])

    #roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])


    #PLOT PREPPING
    fontsize=16
    leg_fontsize=14
    #PLOT 1, LOGMAG
    logmag_plot = plt.figure()
    logmag_ax = plt.subplot(1,1,1)
    plt.ylabel('dB',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    #PLOT 2, Impedance
    impedance_plot = plt.figure()
    impedance_ax1 = plt.subplot(3,1,1)

    plt.ylabel('Re(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax2 = plt.subplot(3,1,2,sharex=impedance_ax1)

    plt.ylabel('Im(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax3 = plt.subplot(3,1,3,sharex=impedance_ax1)

    plt.ylabel('|Z|',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    ignore_roots = ['/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_3.3PF-NO-SHUNT_STRAIGHT_FILLER.csv',\
                    '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF_STRAIGHT_16-75IN.csv',\
                    '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF_STRAIGHT_15-30IN.csv',\
                    '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF_STRAIGHT_135LOGMAG.csv']

    ignore_roots = []

    for root in roots:
        try:
            if root in ignore_roots:
                continue
            print(root)
            if root == '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF-45NH-SHUNT_STRAIGHT_FILLER.csv':
                alpha=1.0
            else:
                alpha=0.8
            logmag_infile = root.replace('_FILLER','_LOGMAG')
            phase_infile = root.replace('_FILLER','_PHASE')

            label = logmag_infile.split('/')[-1].replace('.csv','').replace('_PHASE','').replace('_LOGMAG','').replace('HOURGLASS','HOURGLASS ').replace('_',' ').title().replace('Pf',' pF').replace('Nh',' nH').replace('Ohm',' $\\Omega$').replace('-',' ')
            #cap = logmag_infile.split('/')[-1].split('_')[1]
            #shape = logmag_infile.split('/')[-1].split('_')[2]

            #color = plt.rcParams['axes.prop_cycle'].by_key()['color'][numpy.where(numpy.isin(caps,[cap]))[0][0]]
            #linestyle = linestyles[numpy.where(numpy.isin(shapes,[shape]))[0][0]]

            freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
            freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)
            
            plot_cut_ll = 100            
            plot_cut = freqs/1e6 > plot_cut_ll

            if 'Cat' in label:
                linestyle = '--'
            elif 'Dog' in label:
                linestyle = '-'
            else:
                linestyle = '-.'


            logmag_ax.plot(freqs[plot_cut]/1e6, logmag_vals[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)

            #data = ff.plotLogMag(logmag_infile,freq_range=(min(freqs),max(freqs)),fig=logmag_plot,label=label,color=color,linestyle=linestyle,alpha=0.8)


            lin_vals = ff.logMagToLin(logmag_vals)

            Re, Im = ff.magPhaseToReIm(lin_vals,phase_vals)

            complex_S11 = Re + 1j*Im

            z = ff.linToComplexZ( complex_S11, z_0 = 50.0 )

            impedance_ax1.plot(freqs[plot_cut]/1e6, numpy.real(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax2.plot(freqs[plot_cut]/1e6, numpy.imag(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax3.plot(freqs[plot_cut]/1e6, numpy.abs(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
        except Exception as e:
            print(e)
    logmag_ax.legend(fontsize=leg_fontsize)
    logmag_ax.set_xlim([plot_cut_ll,1000])
    logmag_ax.set_ylim([-30,0.5])


    impedance_ax1.legend(fontsize=leg_fontsize,loc='upper right')


    impedance_ax2.legend(fontsize=leg_fontsize,loc='upper right')

    impedance_ax3.legend(fontsize=leg_fontsize,loc='upper right')

    if mode == 0:
        test_near45(mode=1)

def test():
    '''
    '''
    datapath ='/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/'
    print(datapath)
    plt.rcParams['axes.prop_cycle'].by_key()['color']

    all_linestyles = numpy.array(list(lines.lineStyles.keys()))
    linestyles = numpy.array(['-', '--',':', '-.', ' ', ''], dtype='<U4')#all_linestyles[~numpy.isin(all_linestyles,['None'])]
    
    infiles = numpy.array(glob.glob(datapath + '*100NH*.csv'))
    #infiles = numpy.append(infiles,numpy.array(glob.glob(datapath + '*33*CAT*.csv')))
    
    #infiles = numpy.append(infiles,numpy.array(glob.glob(datapath + '*DOG2*.csv')))
    #infiles = numpy.append(infiles,numpy.array(glob.glob(datapath + '*CAT2*.csv')))
    infiles = numpy.sort(infiles)
    roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])
    print(infiles)

    #caps = numpy.unique([infile.split('/')[-1].split('_')[1] for infile in infiles])
    #shapes = numpy.unique([infile.split('/')[-1].split('_')[2] for infile in infiles])

    #roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])


    #PLOT PREPPING
    fontsize=16
    leg_fontsize=14
    #PLOT 1, LOGMAG
    logmag_plot = plt.figure()
    logmag_ax = plt.subplot(1,1,1)
    plt.ylabel('dB',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    #PLOT 2, Impedance
    impedance_plot = plt.figure()
    impedance_ax1 = plt.subplot(3,1,1)

    plt.ylabel('Re(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax2 = plt.subplot(3,1,2,sharex=impedance_ax1)

    plt.ylabel('Im(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax3 = plt.subplot(3,1,3,sharex=impedance_ax1)

    plt.ylabel('|Z|',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    ignore_roots = ['/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_3.3PF-NO-SHUNT_STRAIGHT_FILLER.csv',\
                    '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF_STRAIGHT_16-75IN.csv',\
                    '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF_STRAIGHT_15-30IN.csv',\
                    '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF_STRAIGHT_135LOGMAG.csv']

    ignore_roots = []

    for root in roots:
        try:
            if root in ignore_roots:
                continue
            print(root)
            if root == '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/HPOL_2PF-45NH-SHUNT_STRAIGHT_FILLER.csv':
                alpha=1.0
            else:
                alpha=0.8
            logmag_infile = root.replace('_FILLER','_LOGMAG')
            phase_infile = root.replace('_FILLER','_PHASE')

            label = logmag_infile.split('/')[-1].replace('.csv','').replace('_PHASE','').replace('_LOGMAG','').replace('HOURGLASS','HOURGLASS ').replace('_',' ').title().replace('Pf',' pF').replace('Nh',' nH').replace('Ohm',' $\\Omega$').replace('-',' ')
            #cap = logmag_infile.split('/')[-1].split('_')[1]
            #shape = logmag_infile.split('/')[-1].split('_')[2]

            #color = plt.rcParams['axes.prop_cycle'].by_key()['color'][numpy.where(numpy.isin(caps,[cap]))[0][0]]
            #linestyle = linestyles[numpy.where(numpy.isin(shapes,[shape]))[0][0]]

            freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
            freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)
            
            plot_cut_ll = 100            
            plot_cut = freqs/1e6 > plot_cut_ll

            if 'Cat' in label:
                linestyle = '--'
            elif 'Dog' in label:
                linestyle = '-'
            else:
                linestyle = '-.'


            logmag_ax.plot(freqs[plot_cut]/1e6, logmag_vals[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)

            #data = ff.plotLogMag(logmag_infile,freq_range=(min(freqs),max(freqs)),fig=logmag_plot,label=label,color=color,linestyle=linestyle,alpha=0.8)


            lin_vals = ff.logMagToLin(logmag_vals)

            Re, Im = ff.magPhaseToReIm(lin_vals,phase_vals)

            complex_S11 = Re + 1j*Im

            z = ff.linToComplexZ( complex_S11, z_0 = 50.0 )

            impedance_ax1.plot(freqs[plot_cut]/1e6, numpy.real(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax2.plot(freqs[plot_cut]/1e6, numpy.imag(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax3.plot(freqs[plot_cut]/1e6, numpy.abs(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
        except Exception as e:
            print(e)
    logmag_ax.legend(fontsize=leg_fontsize)
    logmag_ax.set_xlim([plot_cut_ll,1000])
    logmag_ax.set_ylim([-30,0.5])


    impedance_ax1.legend(fontsize=leg_fontsize,loc='upper right')


    impedance_ax2.legend(fontsize=leg_fontsize,loc='upper right')

    impedance_ax3.legend(fontsize=leg_fontsize,loc='upper right')


def test_2point7():
    '''
    '''
    datapath ='/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_testing_2.7pF/'
    print(datapath)
    plt.rcParams['axes.prop_cycle'].by_key()['color']

    all_linestyles = numpy.array(list(lines.lineStyles.keys()))
    linestyles = numpy.array(['-', '--',':', '-.', ' ', ''], dtype='<U4')#all_linestyles[~numpy.isin(all_linestyles,['None'])]
    
    infiles = numpy.array(glob.glob(datapath + '*.csv'))
    #infiles = numpy.append(infiles,numpy.array(glob.glob(datapath + '*33*CAT*.csv')))

    infiles = numpy.sort(infiles)
    roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])
    print(infiles)

    #caps = numpy.unique([infile.split('/')[-1].split('_')[1] for infile in infiles])
    #shapes = numpy.unique([infile.split('/')[-1].split('_')[2] for infile in infiles])

    #roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])


    #PLOT PREPPING
    fontsize=16
    leg_fontsize=14
    #PLOT 1, LOGMAG
    logmag_plot = plt.figure()
    logmag_ax = plt.subplot(1,1,1)
    plt.ylabel('dB',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    #PLOT 2, Impedance
    impedance_plot = plt.figure()
    impedance_ax1 = plt.subplot(3,1,1)

    plt.ylabel('Re(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax2 = plt.subplot(3,1,2,sharex=impedance_ax1)

    plt.ylabel('Im(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax3 = plt.subplot(3,1,3,sharex=impedance_ax1)

    plt.ylabel('|Z|',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    ignore_roots = ['/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_testing_2.7pF/HPOL_2.7PF-100NH-SHUNT_ALICE-CAT2_FILLER-V2.csv',\
                   '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_testing_2.7pF/HPOL_2.7PF-100NH-SHUNT_ALICE-CAT2_FILLER.csv',\
                   '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_testing_2.7pF/HPOL_2.7PF-100NH-V2_ALICE-CAT2.csv',\
                   '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_testing_2.7pF/HPOL_2.7PF-56NH_ALICE-CAT2_FILLER.csv',\
                   '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_testing_2.7pF/HPOL_2.7PF-NO-SHUNT_ALICE-CAT2_FILLER.csv',\
                   '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_testing_2.7pF/HPOL_2.7PF_47NH-SHUNT_ALICE_CAT2_FILLER.csv',\
                   '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_testing_2.7pF/HPOL_2.7PF-100NH_BOB-DOG2_FILLER.csv']


    for root in roots:
        try:
            if root in ignore_roots:
                continue
            print(root)

            alpha=0.8

            logmag_infile = root.replace('_FILLER','_LOGMAG')
            phase_infile = root.replace('_FILLER','_PHASE')

            label = logmag_infile.split('/')[-1].replace('.csv','').replace('_PHASE','').replace('_LOGMAG','').replace('HOURGLASS','HOURGLASS ').replace('_',' ').title().replace('Pf',' pF').replace('Nh',' nH').replace('Ohm',' $\\Omega$').replace('-',' ')
            
            if 'V3' in label:
                label = label.replace('V3','Missing Capacitor')
            elif 'Shunt' in label:
                label = label.replace('Shunt','Missing Capacitor')
            elif 'V4' in label:
                label = label.replace('V4','All Capacitors')
            elif 'v2' in label:
                label = label.replace('v2','All Capacitors')


            #cap = logmag_infile.split('/')[-1].split('_')[1]
            #shape = logmag_infile.split('/')[-1].split('_')[2]

            #color = plt.rcParams['axes.prop_cycle'].by_key()['color'][numpy.where(numpy.isin(caps,[cap]))[0][0]]
            #linestyle = linestyles[numpy.where(numpy.isin(shapes,[shape]))[0][0]]

            freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
            freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)
            
            plot_cut_ll = 100            
            plot_cut = freqs/1e6 > plot_cut_ll
            '''
            if '100' in label:
                linestyle = '-'
            elif '82' in label:
                linestyle = '-'
            elif '56' in label:
                linestyle = '--'
            elif '47' in label:
                linestyle = '-.'
            else:
                linestyle = ':'

            '''
            linestyle = '-'
            logmag_ax.plot(freqs[plot_cut]/1e6, logmag_vals[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)

            #data = ff.plotLogMag(logmag_infile,freq_range=(min(freqs),max(freqs)),fig=logmag_plot,label=label,color=color,linestyle=linestyle,alpha=0.8)


            lin_vals = ff.logMagToLin(logmag_vals)

            Re, Im = ff.magPhaseToReIm(lin_vals,phase_vals)

            complex_S11 = Re + 1j*Im

            z = ff.linToComplexZ( complex_S11, z_0 = 50.0 )

            impedance_ax1.plot(freqs[plot_cut]/1e6, numpy.real(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax2.plot(freqs[plot_cut]/1e6, numpy.imag(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax3.plot(freqs[plot_cut]/1e6, numpy.abs(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
        except Exception as e:
            print(e)
    logmag_ax.legend(fontsize=leg_fontsize)
    logmag_ax.set_xlim([plot_cut_ll,1000])
    logmag_ax.set_ylim([-30,0.5])


    impedance_ax1.legend(fontsize=leg_fontsize,loc='upper right')


    impedance_ax2.legend(fontsize=leg_fontsize,loc='upper right')

    impedance_ax3.legend(fontsize=leg_fontsize,loc='upper right')









    #PLOT PREPPING
    fontsize=16
    leg_fontsize=14
    #PLOT 1, LOGMAG
    logmag_plot = plt.figure()
    logmag_ax = plt.subplot(1,1,1)
    plt.ylabel('dB',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    #PLOT 2, Impedance
    impedance_plot = plt.figure()
    impedance_ax1 = plt.subplot(3,1,1)

    plt.ylabel('Re(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax2 = plt.subplot(3,1,2,sharex=impedance_ax1)

    plt.ylabel('Im(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax3 = plt.subplot(3,1,3,sharex=impedance_ax1)

    plt.ylabel('|Z|',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    infiles = numpy.array(glob.glob(datapath + '*ALICE*CAT2*.csv'))
    #infiles = numpy.append(infiles,numpy.array(glob.glob(datapath + '*33*CAT*.csv')))

    infiles = numpy.sort(infiles)
    roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])
    print(infiles)



    ignore_roots = ['/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_testing_2.7pF/HPOL_2.7PF-100NH-SHUNT_ALICE-CAT2_FILLER-V2.csv',\
                   '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_testing_2.7pF/HPOL_2.7PF-100NH-SHUNT_ALICE-CAT2_FILLER.csv',\
                   '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_testing_2.7pF/HPOL_2.7PF-100NH-V3_ALICE-CAT2_FILLER.csv',\
                   '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_testing_2.7pF/HPOL_2.7PF-100NH-V2_ALICE-CAT2.csv',\
                   '/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_testing_2.7pF/HPOL_2.7PF-33NH-SHUNT_ALICE-CAT2_FILLER.csv']


    for root in roots:
        try:
            if root in ignore_roots:
                continue
            print(root)

            alpha=0.8

            logmag_infile = root.replace('_FILLER','_LOGMAG')
            phase_infile = root.replace('_FILLER','_PHASE')

            label = logmag_infile.split('/')[-1].replace('.csv','').replace('_PHASE','').replace('_LOGMAG','').replace('HOURGLASS','HOURGLASS ').replace('_',' ').title().replace('Pf',' pF').replace('Nh',' nH').replace('Ohm',' $\\Omega$').replace('-',' ')
            #cap = logmag_infile.split('/')[-1].split('_')[1]
            #shape = logmag_infile.split('/')[-1].split('_')[2]

            #color = plt.rcParams['axes.prop_cycle'].by_key()['color'][numpy.where(numpy.isin(caps,[cap]))[0][0]]
            #linestyle = linestyles[numpy.where(numpy.isin(shapes,[shape]))[0][0]]

            freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
            freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)
            
            plot_cut_ll = 100            
            plot_cut = freqs/1e6 > plot_cut_ll

            '''
            if '100' in label:
                linestyle = '-'
            elif '82' in label:
                linestyle = '-'
            elif '56' in label:
                linestyle = '--'
            elif '47' in label:
                linestyle = '-.'
            else:
                linestyle = ':'

            '''
            linestyle = '-'


            logmag_ax.plot(freqs[plot_cut]/1e6, logmag_vals[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)

            #data = ff.plotLogMag(logmag_infile,freq_range=(min(freqs),max(freqs)),fig=logmag_plot,label=label,color=color,linestyle=linestyle,alpha=0.8)


            lin_vals = ff.logMagToLin(logmag_vals)

            Re, Im = ff.magPhaseToReIm(lin_vals,phase_vals)

            complex_S11 = Re + 1j*Im

            z = ff.linToComplexZ( complex_S11, z_0 = 50.0 )

            impedance_ax1.plot(freqs[plot_cut]/1e6, numpy.real(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax2.plot(freqs[plot_cut]/1e6, numpy.imag(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax3.plot(freqs[plot_cut]/1e6, numpy.abs(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
        except Exception as e:
            print(e)
    logmag_ax.legend(fontsize=leg_fontsize)
    logmag_ax.set_xlim([plot_cut_ll,1000])
    logmag_ax.set_ylim([-30,0.5])


    impedance_ax1.legend(fontsize=leg_fontsize,loc='upper right')


    impedance_ax2.legend(fontsize=leg_fontsize,loc='upper right')

    impedance_ax3.legend(fontsize=leg_fontsize,loc='upper right')

def test_compare():
    '''
    '''
    datapath ='/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_low_resonance/'
    print(datapath)
    plt.rcParams['axes.prop_cycle'].by_key()['color']

    all_linestyles = numpy.array(list(lines.lineStyles.keys()))
    linestyles = numpy.array(['-', '--',':', '-.', ' ', ''], dtype='<U4')#all_linestyles[~numpy.isin(all_linestyles,['None'])]
    
    infiles = numpy.array(glob.glob(datapath + '*'))
    #infiles = numpy.append(infiles,numpy.array(glob.glob(datapath + '*33*CAT*.csv')))

    infiles = numpy.sort(infiles)
    roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])
    print(infiles)

    #caps = numpy.unique([infile.split('/')[-1].split('_')[1] for infile in infiles])
    #shapes = numpy.unique([infile.split('/')[-1].split('_')[2] for infile in infiles])

    #roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])


    #PLOT PREPPING
    fontsize=16
    leg_fontsize=14
    #PLOT 1, LOGMAG
    logmag_plot = plt.figure()
    logmag_ax = plt.subplot(1,1,1)
    plt.ylabel('dB',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    #PLOT 2, Impedance
    impedance_plot = plt.figure()
    impedance_ax1 = plt.subplot(3,1,1)

    plt.ylabel('Re(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax2 = plt.subplot(3,1,2,sharex=impedance_ax1)

    plt.ylabel('Im(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax3 = plt.subplot(3,1,3,sharex=impedance_ax1)

    plt.ylabel('|Z|',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    ignore_roots = []


    for root in roots:
        try:
            if root in ignore_roots:
                continue
            print(root)

            alpha=0.8

            if '_CP-' in root:
                logmag_infile = root.replace('_FILLER','_LOGMAG')
                #phase_infile = root.replace('_FILLER','_PHASE')

                label = logmag_infile.split('/')[-1].replace('.csv','').replace('.txt','').replace('_PHASE','').replace('_LOGMAG','').replace('HOURGLASS','HOURGLASS ').replace('_',' ').title().replace('Pf',' pF').replace('Nh',' nH').replace('Ohm',' $\\Omega$').replace('-',' ')
                
                if 'Cp' in label:
                    label = label.replace('Cp','\nCal Poly')


                #cap = logmag_infile.split('/')[-1].split('_')[1]
                #shape = logmag_infile.split('/')[-1].split('_')[2]

                #color = plt.rcParams['axes.prop_cycle'].by_key()['color'][numpy.where(numpy.isin(caps,[cap]))[0][0]]
                #linestyle = linestyles[numpy.where(numpy.isin(shapes,[shape]))[0][0]]

                freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=8)
                phase_vals = numpy.zeros_like(logmag_vals)
            else:

                logmag_infile = root.replace('_FILLER','_LOGMAG')
                phase_infile = root.replace('_FILLER','_PHASE')

                label = logmag_infile.split('/')[-1].replace('.csv','').replace('_PHASE','').replace('_LOGMAG','').replace('HOURGLASS','HOURGLASS ').replace('_',' ').title().replace('Pf',' pF').replace('Nh',' nH').replace('Ohm',' $\\Omega$').replace('-',' ')
                
                if 'V3' in label:
                    label = label.replace('V3','Missing Capacitor')
                elif 'Shunt' in label:
                    label = label.replace('Shunt','Missing Capacitor')
                elif 'V4' in label:
                    label = label.replace('V4','All Capacitors')
                elif 'v2' in label:
                    label = label.replace('v2','All Capacitors')


                #cap = logmag_infile.split('/')[-1].split('_')[1]
                #shape = logmag_infile.split('/')[-1].split('_')[2]

                #color = plt.rcParams['axes.prop_cycle'].by_key()['color'][numpy.where(numpy.isin(caps,[cap]))[0][0]]
                #linestyle = linestyles[numpy.where(numpy.isin(shapes,[shape]))[0][0]]

                freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
                freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)



            
            plot_cut_ll = 100            
            plot_cut = freqs/1e6 > plot_cut_ll

            linestyle = '-'


            logmag_ax.plot(freqs[plot_cut]/1e6, logmag_vals[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)

            #data = ff.plotLogMag(logmag_infile,freq_range=(min(freqs),max(freqs)),fig=logmag_plot,label=label,color=color,linestyle=linestyle,alpha=0.8)


            lin_vals = ff.logMagToLin(logmag_vals)

            Re, Im = ff.magPhaseToReIm(lin_vals,phase_vals)

            complex_S11 = Re + 1j*Im

            z = ff.linToComplexZ( complex_S11, z_0 = 50.0 )

            impedance_ax1.plot(freqs[plot_cut]/1e6, numpy.real(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax2.plot(freqs[plot_cut]/1e6, numpy.imag(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax3.plot(freqs[plot_cut]/1e6, numpy.abs(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
        except Exception as e:
            print(e)
    logmag_ax.legend(fontsize=leg_fontsize)
    logmag_ax.set_xlim([plot_cut_ll,1000])
    logmag_ax.set_ylim([-30,0.5])


    impedance_ax1.legend(fontsize=leg_fontsize,loc='upper right')


    impedance_ax2.legend(fontsize=leg_fontsize,loc='upper right')

    impedance_ax3.legend(fontsize=leg_fontsize,loc='upper right')


def test_caps():
    '''
    '''
    datapath ='/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/compare_caps/'
    print(datapath)
    plt.rcParams['axes.prop_cycle'].by_key()['color']

    all_linestyles = numpy.array(list(lines.lineStyles.keys()))
    linestyles = numpy.array(['-', '--',':', '-.', ' ', ''], dtype='<U4')#all_linestyles[~numpy.isin(all_linestyles,['None'])]
    
    infiles = numpy.array(glob.glob(datapath + '*'))
    #infiles = numpy.append(infiles,numpy.array(glob.glob(datapath + '*33*CAT*.csv')))

    infiles = numpy.sort(infiles)
    roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])
    print(infiles)

    #caps = numpy.unique([infile.split('/')[-1].split('_')[1] for infile in infiles])
    #shapes = numpy.unique([infile.split('/')[-1].split('_')[2] for infile in infiles])

    #roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])


    #PLOT PREPPING
    fontsize=16
    leg_fontsize=14
    #PLOT 1, LOGMAG
    logmag_plot = plt.figure()
    logmag_ax = plt.subplot(1,1,1)
    plt.ylabel('dB',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    #PLOT 2, Impedance
    impedance_plot = plt.figure()
    impedance_ax1 = plt.subplot(3,1,1)

    plt.ylabel('Re(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax2 = plt.subplot(3,1,2,sharex=impedance_ax1)

    plt.ylabel('Im(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax3 = plt.subplot(3,1,3,sharex=impedance_ax1)

    plt.ylabel('|Z|',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    ignore_roots = []


    for root in roots:
        try:
            if root in ignore_roots:
                continue
            print(root)

            alpha=0.8

            if '_CP-' in root:
                logmag_infile = root.replace('_FILLER','_LOGMAG')
                #phase_infile = root.replace('_FILLER','_PHASE')

                label = logmag_infile.split('/')[-1].replace('.csv','').replace('.txt','').replace('_PHASE','').replace('_LOGMAG','').replace('HOURGLASS','HOURGLASS ').replace('_',' ').title().replace('Pf',' pF').replace('Nh',' nH').replace('Ohm',' $\\Omega$').replace('-',' ')
                
                if 'Cp' in label:
                    label = label.replace('Cp','\nCal Poly')


                #cap = logmag_infile.split('/')[-1].split('_')[1]
                #shape = logmag_infile.split('/')[-1].split('_')[2]

                #color = plt.rcParams['axes.prop_cycle'].by_key()['color'][numpy.where(numpy.isin(caps,[cap]))[0][0]]
                #linestyle = linestyles[numpy.where(numpy.isin(shapes,[shape]))[0][0]]

                freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=8)
                phase_vals = numpy.zeros_like(logmag_vals)
            else:

                logmag_infile = root.replace('_FILLER','_LOGMAG')
                phase_infile = root.replace('_FILLER','_PHASE')

                label = logmag_infile.split('/')[-1].replace('.csv','').replace('_PHASE','').replace('_LOGMAG','').replace('HOURGLASS','HOURGLASS ').replace('_',' ').title().replace('Pf',' pF').replace('Nh',' nH').replace('Ohm',' $\\Omega$').replace('-',' ')
                
                if 'V3' in label:
                    label = label.replace('V3','Missing Capacitor')
                elif 'Shunt' in label:
                    label = label.replace('Shunt','Missing Capacitor')
                elif 'V4' in label:
                    label = label.replace('V4','All Capacitors')
                elif 'v2' in label:
                    label = label.replace('v2','All Capacitors')


                #cap = logmag_infile.split('/')[-1].split('_')[1]
                #shape = logmag_infile.split('/')[-1].split('_')[2]

                #color = plt.rcParams['axes.prop_cycle'].by_key()['color'][numpy.where(numpy.isin(caps,[cap]))[0][0]]
                #linestyle = linestyles[numpy.where(numpy.isin(shapes,[shape]))[0][0]]

                freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
                freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)



            
            plot_cut_ll = 100            
            plot_cut = freqs/1e6 > plot_cut_ll

            linestyle = '-'


            logmag_ax.plot(freqs[plot_cut]/1e6, logmag_vals[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)

            #data = ff.plotLogMag(logmag_infile,freq_range=(min(freqs),max(freqs)),fig=logmag_plot,label=label,color=color,linestyle=linestyle,alpha=0.8)


            lin_vals = ff.logMagToLin(logmag_vals)

            Re, Im = ff.magPhaseToReIm(lin_vals,phase_vals)

            complex_S11 = Re + 1j*Im

            z = ff.linToComplexZ( complex_S11, z_0 = 50.0 )

            impedance_ax1.plot(freqs[plot_cut]/1e6, numpy.real(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax2.plot(freqs[plot_cut]/1e6, numpy.imag(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
            impedance_ax3.plot(freqs[plot_cut]/1e6, numpy.abs(z)[plot_cut],label=label,alpha=alpha,linestyle=linestyle)#,color=color)
        except Exception as e:
            print(e)
    logmag_ax.legend(fontsize=leg_fontsize)
    logmag_ax.set_xlim([plot_cut_ll,1000])
    logmag_ax.set_ylim([-30,0.5])


    impedance_ax1.legend(fontsize=leg_fontsize,loc='upper right')


    impedance_ax2.legend(fontsize=leg_fontsize,loc='upper right')

    impedance_ax3.legend(fontsize=leg_fontsize,loc='upper right')

if __name__ == '__main__':

    test_caps()


    #test_near45()
    #test_switch_to_consistent_wire_lengths()
    # test_2point7()
    # test_compare()
    '''

    ################################
    #COMPARING SHIELDED AND UNSHIELDED
    ################################

    datapath ='/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_data/shield/'
    infiles = glob.glob(datapath + '*.csv')

    roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])


    #PLOT PREPPING

    #PLOT 1, LOGMAG
    logmag_plot = plt.figure()
    logmag_ax = plt.subplot(1,1,1)
    plt.ylabel('dB',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)


    #PLOT 2, Impedance
    impedance_plot = plt.figure()
    impedance_ax1 = plt.subplot(3,1,1)

    plt.ylabel('Re(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax2 = plt.subplot(3,1,2)

    plt.ylabel('Im(Z)',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    impedance_ax3 = plt.subplot(3,1,3)

    plt.ylabel('|Z|',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    for root in roots:
        logmag_infile = root.replace('_FILLER','_LOGMAG')
        phase_infile = root.replace('_FILLER','_PHASE')

        label = logmag_infile.split('/')[-1].replace('.csv','').replace('_PHASE','').replace('_LOGMAG','').replace('HOURGLASS','HOURGLASS ').replace('_',' ').title().replace('Pf',' pF').replace('Nh',' nH').replace('Ohm',' $\\Omega$')
        cap = logmag_infile.split('/')[-1].split('_')[1]
        shape = logmag_infile.split('/')[-1].split('_')[2]

        color = plt.rcParams['axes.prop_cycle'].by_key()['color'][numpy.where(numpy.isin(caps,[cap]))[0][0]]
        linestyle = linestyles[numpy.where(numpy.isin(shapes,[shape]))[0][0]]

        freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
        freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)
        
        plot_cut = freqs/1e6 > 200


        logmag_ax.plot(freqs[plot_cut]/1e6, logmag_vals[plot_cut],label=label,linestyle=linestyle,color=color)

        #data = ff.plotLogMag(logmag_infile,freq_range=(min(freqs),max(freqs)),fig=logmag_plot,label=label,color=color,linestyle=linestyle,alpha=0.8)


        lin_vals = ff.logMagToLin(logmag_vals)

        Re, Im = ff.magPhaseToReIm(lin_vals,phase_vals)

        complex_S11 = Re + 1j*Im

        z = ff.linToComplexZ( complex_S11, z_0 = 50.0 )

        impedance_ax1.plot(freqs[plot_cut]/1e6, numpy.real(z)[plot_cut],label=label,linestyle=linestyle,color=color)
        impedance_ax2.plot(freqs[plot_cut]/1e6, numpy.imag(z)[plot_cut],label=label,linestyle=linestyle,color=color)
        impedance_ax3.plot(freqs[plot_cut]/1e6, numpy.abs(z)[plot_cut],label=label,linestyle=linestyle,color=color)





    logmag_ax.legend(fontsize=leg_fontsize)
    logmag_ax.set_xlim([200,1000])


    impedance_ax1.legend(fontsize=leg_fontsize,loc='upper right')


    impedance_ax2.legend(fontsize=leg_fontsize,loc='upper right')

    impedance_ax3.legend(fontsize=leg_fontsize,loc='upper right')

    '''