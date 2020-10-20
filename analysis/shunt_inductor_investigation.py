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

def makeGeneralPlots2():
    '''
    This will make s11 expected v.s. measured plots in both logmag and smith chart.   As well as the general
    curve for expected points below threshold.  
    '''
    try:
        datapath ='/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_copper_tab_testing_aug2020/s11/'
        starter_root_match = 'hpol_ape_110nh_2p7pf' #NOSHUNT
        starter_cap_value = 2.7e-12 #farads
        starter_shunt_inductor_value = 110.0e-9 #Henries
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
        feed = SmithMatcher(starter_freqs, starter_logmag_vals, starter_phase_vals, initial_capacitor_value=starter_cap_value, initial_shunt_inductor_value=starter_shunt_inductor_value, z_0=50)

        #inductor_values_nH = numpy.linspace(15,250,300)
        inductor_values_of_interest = numpy.array([])
        inductor_values_nH = numpy.linspace(5.0,200.0,500)#numpy.array([47.0,56,100,110])
        test_statistic = numpy.zeros_like(inductor_values_nH)
        db_cut = -4 #dB


        for index, L in enumerate(inductor_values_nH):
            feed.swapRLC('pl',L*1.0e-9)
            test_statistic[index] = feed.getPercentBelow(db=db_cut)
            for roi in inductor_values_of_interest:
                if numpy.argmin(numpy.abs(roi - inductor_values_nH)) == index:
                    feed.plotCurrentLogMagS11(label='L = %0.2f, (closest to requested %0.2f)'%(L,roi))
            feed.reset()


        fig = plt.figure()
        curve_ax = plt.gca()
        plt.plot(inductor_values_nH, test_statistic,label='Predicted Curve from Varying %.2f nH Feed'%(starter_shunt_inductor_value*1e9))
        plt.xlabel('Individual Added Shunt Inductor Values (nH)')
        plt.ylabel('Fraction of S11 Points below %0.2f dB'%db_cut)
        plt.minorticks_on()
        plt.grid(b=True, which='major', color='k', linestyle='-',alpha=0.75)
        plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.25)
        feed.reset()
        plt.axhline(feed.getPercentBelow(db=db_cut),c='r',linestyle='--',label='Measured %.2f nH Fraction'%(starter_shunt_inductor_value*1e9))

        measured_test_statistic = numpy.zeros(len(roots))
        measured_inductor_values_nH = numpy.zeros(len(roots))
        for index, root in enumerate(roots):
            logmag_infile = root.replace('_FILLER','_LOGMAG')
            phase_infile = root.replace('_FILLER','_PHASE')
            L = root.split('/')[-1].split('_')[2] #nH
            if L.lower() == 'noshunt':
                L = 0
            else:
                L = float(L.lower().replace('nh',''))
            C = float(root.split('/')[-1].split('_')[3].lower().replace('pf','').replace('p','.')) #pF
            measured_inductor_values_nH[index] = L
            freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
            freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)
            measured_feed = SmithMatcher(freqs, logmag_vals, phase_vals, initial_capacitor_value=C*1e-12, initial_shunt_inductor_value=L*1e-9, z_0=50)
            measured_test_statistic[index] = measured_feed.getPercentBelow(db=db_cut)
            #measured_feed.plotCurrentLogMagS11(label=root)
            
            feed.reset()
            compare_fig = plt.figure()
            compare_ax1 = plt.subplot(2,1,1)

            compare_ax1 = feed.plotCurrentLogMagS11(ax=compare_ax1,label='Measured: %0.2f nH Inductors'%(starter_shunt_inductor_value*1e9)) #Labelling assumes starter_ curves are for no shunt curve.
            feed.swapRLC('pl',L*1e-9)
            compare_ax1 = measured_feed.plotCurrentLogMagS11(ax=compare_ax1, label='Measured: %0.2f nH Inductors'%L)
            compare_ax1 = feed.plotCurrentLogMagS11(ax=compare_ax1, label='Predicted: %0.2f nH Inductors'%L)
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
            compare_ax2 = feed.plotSmithChart(ax=compare_ax2) #Labelling assumes starter_ curves are for no shunt curve.
            feed.swapRLC('pl',L*1e-9)
            compare_ax2 = measured_feed.plotSmithChart(ax=compare_ax2)
            compare_ax2 = feed.plotSmithChart(ax=compare_ax2)



        curve_ax.scatter(measured_inductor_values_nH,measured_test_statistic,color='r',edgecolors='k',label='Measured Feed Values')
        curve_ax.legend(loc='upper right')

    except Exception as e:
        print('\nError in %s'%inspect.stack()[0][3])
        print(e)
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)

def makeGeneralPlots3():
    '''
    This will make s11 expected v.s. measured plots in both logmag and smith chart.   As well as the general
    curve for expected points below threshold.  
    '''
    readout_transformer_factor = 1.0
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
        inductor_values_of_interest = numpy.array([])
        inductor_values_nH = numpy.linspace(1.0,100.0,50)#numpy.array([47.0,56,100,110])
        test_statistic = numpy.zeros_like(inductor_values_nH)
        
        cap_values_pF = numpy.linspace(1.0,10.0,50)
        cap_values_of_interest = numpy.array([])
        test_statistic_cap = numpy.zeros_like(cap_values_pF)
        fontsize=16


        db_cut = -4 #dB
        freq_low_cut_MHz=200
        freq_high_cut_MHz=800
        test_C = 4.7 #pF
        test_L = 47 #nH



        #Prep calculations
        starter_shunt_inductor_value = starter_root_match.split('/')[-1].split('_')[2] #nH
        if starter_shunt_inductor_value.lower() == 'noshunt':
            starter_shunt_inductor_value = 0
        else:
            starter_shunt_inductor_value = float(starter_shunt_inductor_value.lower().replace('nh',''))*1e-9 #H
        starter_cap_value = float(starter_root_match.split('/')[-1].split('_')[3].lower().replace('pf','').replace('p','.'))*1e-12 #F
            
        mesh_L, mesh_C = numpy.meshgrid(inductor_values_nH,cap_values_pF)
        mesh_ts = numpy.zeros_like(mesh_L)

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


        for c_index, C in enumerate(cap_values_pF):
            sys.stdout.write('(%i/%i)\t\t\t\r'%(c_index+1,len(cap_values_pF)))
            sys.stdout.flush()
            for l_index, L in enumerate(inductor_values_nH):
                feed.swapRLC('sc',C*1.0e-12)
                feed.swapRLC('pl',L*1.0e-9)
                mesh_ts[c_index][l_index] = feed.getPercentBelow(db=db_cut,freq_low_cut_MHz=freq_low_cut_MHz,freq_high_cut_MHz=freq_high_cut_MHz)
                feed.reset()

        fig = plt.figure()
        mesh_ax = plt.gca()
        im = mesh_ax.pcolormesh(mesh_L, mesh_C, mesh_ts, vmin=None, vmax=None,cmap=plt.cm.coolwarm)
        cbar = fig.colorbar(im)
        cbar.set_label('Fraction of S11 Points below %0.2f dB [%i - %i MHz]'%(db_cut,freq_low_cut_MHz,freq_high_cut_MHz),fontsize=fontsize)
        plt.xlabel('Individual Swapped Shunt Inductor Values (nH)',fontsize=fontsize)
        plt.ylabel('Individual Swapped Series Cap Values (pF)',fontsize=fontsize)
        #plt.zlabel('Fraction of S11 Points below %0.2f dB'%db_cut)
        plt.minorticks_on()
        plt.grid(b=True, which='major', color='k', linestyle='-',alpha=0.75)
        plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.25)
        plt.axhline(2.7,linestyle='--',c='k')
        feed.reset()

        #TEMPORARY TEST OF CAP VALUES 
        feed.reset()
        feed.swapRLC('sc',test_C*1.0e-12)
        feed.swapRLC('pl',test_L*1.0e-9)
        feed.plotCurrentLogMagS11(label='C = %0.2f, L = %0.2f '%(test_C,test_L),plot_cut_ll=0,plot_cut_ul=1500)
        feed.plotSmithChart(label='C = %0.2f, L = %0.2f '%(test_C,test_L),plot_cut_ll=0,plot_cut_ul=1500)
        feed.reset()

    
        for index, C in enumerate(cap_values_pF):
            feed.swapRLC('sc',C*1.0e-12)
            test_statistic_cap[index] = feed.getPercentBelow(db=db_cut,freq_low_cut_MHz=freq_low_cut_MHz,freq_high_cut_MHz=freq_high_cut_MHz)
            for roi in cap_values_of_interest:
                if numpy.argmin(numpy.abs(roi - cap_values_pF)) == index:
                    feed.plotCurrentLogMagS11(label='C = %0.2f, (closest to requested %0.2f)'%(C,roi))
            feed.reset()

        for index, L in enumerate(inductor_values_nH):
            feed.swapRLC('pl',L*1.0e-9)
            test_statistic[index] = feed.getPercentBelow(db=db_cut,freq_low_cut_MHz=freq_low_cut_MHz,freq_high_cut_MHz=freq_high_cut_MHz)
            for roi in inductor_values_of_interest:
                if numpy.argmin(numpy.abs(roi - inductor_values_nH)) == index:
                    feed.plotCurrentLogMagS11(label='L = %0.2f, (closest to requested %0.2f)'%(L,roi))
            feed.reset()


        fig = plt.figure()
        curve_ax = plt.gca()
        plt.plot(cap_values_pF, test_statistic_cap,label='Predicted Curve from Varying %.2f pF %.2f nH Feed'%(starter_cap_value*1e12,starter_shunt_inductor_value*1e9))
        plt.xlabel('Individual Swapped Series Cap Values (pF)')
        plt.ylabel('Fraction of S11 Points below %0.2f dB'%db_cut)
        plt.minorticks_on()
        plt.grid(b=True, which='major', color='k', linestyle='-',alpha=0.75)
        plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.25)
        feed.reset()


        fig = plt.figure()
        curve_ax = plt.gca()
        plt.plot(inductor_values_nH, test_statistic,label='Predicted Curve from Varying %.2f nH Feed'%(starter_shunt_inductor_value*1e9))
        plt.xlabel('Individual Added Shunt Inductor Values (nH)')
        plt.ylabel('Fraction of S11 Points below %0.2f dB [%i - %i MHz]'%(db_cut,freq_low_cut_MHz,freq_high_cut_MHz))
        plt.minorticks_on()
        plt.grid(b=True, which='major', color='k', linestyle='-',alpha=0.75)
        plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.25)
        feed.reset()
        plt.axhline(feed.getPercentBelow(db=db_cut,freq_low_cut_MHz=freq_low_cut_MHz,freq_high_cut_MHz=freq_high_cut_MHz),c='r',linestyle='--',label='Measured %.2f nH Fraction'%(starter_shunt_inductor_value*1e9))

        measured_test_statistic = numpy.zeros(len(roots))
        measured_inductor_values_nH = numpy.zeros(len(roots))
        for index, root in enumerate(roots):
            logmag_infile = root.replace('_FILLER','_LOGMAG')
            phase_infile = root.replace('_FILLER','_PHASE')
            L = root.split('/')[-1].split('_')[2] #nH
            if L.lower() == 'noshunt':
                L = 0
            else:
                L = float(L.lower().replace('nh',''))
            C = float(root.split('/')[-1].split('_')[3].lower().replace('pf','').replace('p','.')) #pF
            measured_inductor_values_nH[index] = L
            freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
            freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)
            measured_feed = SmithMatcher(freqs, logmag_vals, phase_vals, initial_capacitor_value=C*1e-12, initial_shunt_inductor_value=L*1e-9, z_0=50, readout_transformer_factor=readout_transformer_factor)
            measured_test_statistic[index] = measured_feed.getPercentBelow(db=db_cut,freq_low_cut_MHz=freq_low_cut_MHz,freq_high_cut_MHz=freq_high_cut_MHz)
            #measured_feed.plotCurrentLogMagS11(label=root)
            
            feed.reset()
            compare_fig = plt.figure()
            compare_ax1 = plt.subplot(2,1,1)

            compare_ax1 = feed.plotCurrentLogMagS11(ax=compare_ax1,label='Measured: %0.2f nH Inductors'%(starter_shunt_inductor_value*1e9)) #Labelling assumes starter_ curves are for no shunt curve.
            feed.swapRLC('pl',L*1e-9)
            compare_ax1 = measured_feed.plotCurrentLogMagS11(ax=compare_ax1, label='Measured: %0.2f nH Inductors'%L)
            compare_ax1 = feed.plotCurrentLogMagS11(ax=compare_ax1, label='Predicted: %0.2f nH Inductors'%L)
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
            compare_ax2 = feed.plotSmithChart(ax=compare_ax2) #Labelling assumes starter_ curves are for no shunt curve.
            feed.swapRLC('pl',L*1e-9)
            compare_ax2 = measured_feed.plotSmithChart(ax=compare_ax2)
            compare_ax2 = feed.plotSmithChart(ax=compare_ax2)



        curve_ax.scatter(measured_inductor_values_nH,measured_test_statistic,color='r',edgecolors='k',label='Measured Feed Values')
        curve_ax.legend(loc='upper right')

    except Exception as e:
        print('\nError in %s'%inspect.stack()[0][3])
        print(e)
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)

def makeGeneralPlots():
    '''
    This will make s11 expected v.s. measured plots in both logmag and smith chart.   As well as the general
    curve for expected points below threshold.  
    '''
    try:
        #datapath ='/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_copper_tab_testing_aug2020/s11/'
        datapath ='/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_al_tube_testing_sep2020/s11/'
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
        inductor_values_of_interest = numpy.array([])
        inductor_values_nH = numpy.linspace(5.0,200.0,500)#numpy.array([47.0,56,100,110])
        test_statistic = numpy.zeros_like(inductor_values_nH)
        db_cut = -4 #dB

        #L wiLL given as an "Effective Inductance" given by L/effective_inductor_factor. 
        effective_inductor_factor = 3.0 #We thought this should be 3 but testing with 1 makes things more accurate to measurements.  

        for index, L in enumerate(inductor_values_nH):
            feed.addRLC('pl',(L/effective_inductor_factor)*1.0e-9)
            test_statistic[index] = feed.getPercentBelow(db=db_cut)
            for roi in inductor_values_of_interest:
                if numpy.argmin(numpy.abs(roi - inductor_values_nH)) == index:
                    feed.plotCurrentLogMagS11(label='L = %0.2f, (closest to requested %0.2f)'%(L,roi))
            feed.reset()


        fig = plt.figure()
        curve_ax = plt.gca()
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
            
            feed.reset()
            compare_fig = plt.figure()
            compare_ax1 = plt.subplot(2,1,1)

            compare_ax1 = feed.plotCurrentLogMagS11(ax=compare_ax1,label='Measured: No Inductors') #Labelling assumes starter_ curves are for no shunt curve.
            feed.addRLC('pl',(L/effective_inductor_factor)*1e-9)
            compare_ax1 = measured_feed.plotCurrentLogMagS11(ax=compare_ax1, label='Measured: %0.2f nH Inductors'%L)
            compare_ax1 = feed.plotCurrentLogMagS11(ax=compare_ax1, label='Predicted: %0.2f nH Inductors'%L)
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
            compare_ax2 = feed.plotSmithChart(ax=compare_ax2) #Labelling assumes starter_ curves are for no shunt curve.
            feed.addRLC('pl',(L/effective_inductor_factor)*1e-9)
            compare_ax2 = measured_feed.plotSmithChart(ax=compare_ax2)
            compare_ax2 = feed.plotSmithChart(ax=compare_ax2)



        curve_ax.scatter(measured_inductor_values_nH,measured_test_statistic,color='r',edgecolors='k',label='Measured Feed Values')
        curve_ax.legend(loc='upper right')

    except Exception as e:
        rint('\nError in %s'%inspect.stack()[0][3])
        print(e)
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)


def testing():
    try:
        #

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

        #L wiLL given as an "Effective Inductance" given by L/effective_inductor_factor. 
        #effective_inductor_factor = 1.0 #We thought this should be 3 but testing with 1 makes things more accurate to measurements.  
        measured_inductor_values_nH = numpy.zeros(len(roots))

        plt.figure()
        w = numpy.ones(len(feed.freqs))
        w[feed.freqs/1e6 < 450] = 3
        w[numpy.logical_and(feed.freqs/1e6 >= 450,feed.freqs/1e6 < 550)] = numpy.linspace(1,3,sum(numpy.logical_and(feed.freqs/1e6 >= 450,feed.freqs/1e6 < 550)))[::-1]
        plt.plot(feed.freqs,w)
        plt.plot(feed.freqs,4 - numpy.linspace(1,numpy.sqrt(3),len(feed.freqs))**2)
        plt.plot(feed.freqs,4 - numpy.linspace(1,3**(1/0.25),len(feed.freqs))**0.25)
        plt.plot(feed.freqs,4 - numpy.linspace(1,3,len(feed.freqs)))

        test_resistor_values = [0.0]#[0.0,0.5/3,0.5,1.5]
        #for effective_inductor_factor in [1.0,3.0,4 - numpy.linspace(1,numpy.sqrt(3),len(feed.freqs))**2]:
        for effective_inductor_factor in [w]:
            for plot_cut_ll, plot_cut_ul in [[250,500],[500,850],[100,850]]:
                for index, root in enumerate(roots):
                    logmag_infile = root.replace('_FILLER','_LOGMAG')
                    phase_infile = root.replace('_FILLER','_PHASE')
                    L = root.split('/')[-1].split('_')[2]
                    if L.lower() == 'noshunt':
                        L = 0
                    else:
                        L = float(L.lower().replace('nh',''))
                    if L != 22:
                        continue
                    measured_inductor_values_nH[index] = L
                    freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
                    freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)
                    measured_feed = Feed(freqs, logmag_vals, phase_vals, z_0=50)
                    #measured_feed.plotCurrentLogMagS11(label=root)
                    
                    feed.reset()
                    compare_fig = plt.figure()
                    plt.suptitle('%0.1f to %0.1f MHz'%(plot_cut_ll,plot_cut_ul))
                    compare_ax1 = plt.subplot(2,1,1)

                    compare_ax1 = feed.plotCurrentLogMagS11(ax=compare_ax1,label='Measured: No Inductors',plot_cut_ll=plot_cut_ll,plot_cut_ul=plot_cut_ul) #Labelling assumes starter_ curves are for no shunt curve.
                    for added_series_resistance in test_resistor_values:
                        if added_series_resistance != 0:
                            feed.addRLC('sr',added_series_resistance)
                        #feed.addRLC('pl',(L/effective_inductor_factor)*1e-9)
                        feed.add22nHShuntInductor(effective_inductor_factor=effective_inductor_factor)
                        compare_ax1 = feed.plotCurrentLogMagS11(ax=compare_ax1, label='Predicted: %0.2f nH Inductors\nAdded %0.4f Ohm Series R'%(L,added_series_resistance),plot_cut_ll=plot_cut_ll,plot_cut_ul=plot_cut_ul)
                        feed.reset()

                    compare_ax1 = measured_feed.plotCurrentLogMagS11(ax=compare_ax1, label='Measured: %0.2f nH Inductors'%L,plot_cut_ll=plot_cut_ll,plot_cut_ul=plot_cut_ul)
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
                    compare_ax2 = feed.plotSmithChart(ax=compare_ax2,plot_cut_ll=plot_cut_ll,plot_cut_ul=plot_cut_ul) #Labelling assumes starter_ curves are for no shunt curve.
                    for added_series_resistance in test_resistor_values:
                        if added_series_resistance != 0:
                            feed.addRLC('sr',added_series_resistance)                    
                        #feed.addRLC('pl',(L/effective_inductor_factor)*1e-9)
                        feed.add22nHShuntInductor(effective_inductor_factor=effective_inductor_factor)
                        compare_ax2 = feed.plotSmithChart(ax=compare_ax2,plot_cut_ll=plot_cut_ll,plot_cut_ul=plot_cut_ul)
                        feed.reset()
                    compare_ax2 = measured_feed.plotSmithChart(ax=compare_ax2,plot_cut_ll=plot_cut_ll,plot_cut_ul=plot_cut_ul)


    except Exception as e:
        print(e)
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)

def matchFeed1toFeed2(feed1,feed2,expected_L_range=[5,150],n=100,plot=False):
    '''
    This will predict the inductor value needed at each frequency to move the points from feed 1 onto feed 2.
    Will use the value of feed1 and feed2 that occur when they are reset.
    '''
    try:
        feed1.reset()
        feed2.reset()
        inductor_attempted_values = numpy.linspace(min(expected_L_range),max(expected_L_range),n)
        output_inductor_values = numpy.zeros_like(feed1.freqs)
        current_best_match = numpy.ones_like(feed1.freqs)*1e9 #These are the absolute difference between the adjusted feed1 and feed2.  Used for comparing if new value is best.
        for L in inductor_attempted_values:
            feed1.addRLC('pl',L*1e-9)
            diff = numpy.abs((numpy.angle(feed1.adjusted_complexs11) - numpy.angle(feed2.adjusted_complexs11)%(2*numpy.pi)))
            improved_cut = diff < current_best_match
            current_best_match[improved_cut] = diff[improved_cut]
            output_inductor_values[improved_cut] = L
            feed1.reset()
        if plot == True:
            plt.figure()
            plt.plot(feed1.freqs/1e6, output_inductor_values)
            plt.ylabel('Inductor Value for Best Match')
            plt.xlabel('Freqs (MHz)')
        return feed1.freqs, output_inductor_values
    except Exception as e:
        print('\nError in %s'%inspect.stack()[0][3])
        print(e)
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)



if __name__ == '__main__':
    #plt.close('all')
    #makeGeneralPlots2()
    makeGeneralPlots3()
    
    '''try:
                    #
            
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
                    feed1 = Feed(starter_freqs, starter_logmag_vals, starter_phase_vals, z_0=50)
            
                    #Prepare the starting S11 that is to be shifted around.
                    for root in roots:
                        L = root.split('/')[-1].split('_')[2]
                        if L.lower() == 'noshunt':
                            L = 0
                        else:
                            L = float(L.lower().replace('nh',''))
                        if L == 22: 
                            break
            
                    logmag_infile = root.replace('_FILLER','_LOGMAG')
                    phase_infile = root.replace('_FILLER','_PHASE')
            
                    #Load intial data and make feed
                    freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
                    freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)
                    feed2 = Feed(freqs, logmag_vals, phase_vals, z_0=50)
            
                    out_freqs, out_L = matchFeed1toFeed2(feed1,feed2,expected_L_range=[22/3,22],n=100,plot=True)
            
            
                    feed1.reset()
                    compare_fig = plt.figure()
                    compare_ax1 = plt.subplot(2,1,1)
            
                    compare_ax1 = feed1.plotCurrentLogMagS11(ax=compare_ax1,label='Measured: No Inductors') #Labelling assumes starter_ curves are for no shunt curve.
                    feed1.add22nHShuntInductor(effective_inductor_factor=22/out_L)
                    compare_ax1 = feed2.plotCurrentLogMagS11(ax=compare_ax1, label='Measured: %0.2f nH Inductors'%L)
                    compare_ax1 = feed1.plotCurrentLogMagS11(ax=compare_ax1, label='Variable Inductor near 22 nH')
                    plt.ylabel('dB',fontsize=14)
                    plt.xlabel('MHz',fontsize=14)
                    plt.minorticks_on()
                    plt.grid(b=True, which='major', color='k', linestyle='-')
                    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)
            
                    feed1.reset()
                    #compare_ax2 = plt.subplot(2,1,2)
                    compare_ax2 = get_smith(compare_fig, rect = 212, plot_impedance = False, plot_ticks = True, plot_admittance = True, plot_labels = False)
                    plt.ylabel('Im($\\Gamma$)',fontsize=14)
                    plt.xlabel('Re($\\Gamma$)',fontsize=14)
                    compare_ax2 = feed1.plotSmithChart(ax=compare_ax2) #Labelling assumes starter_ curves are for no shunt curve.
                    feed1.add22nHShuntInductor(effective_inductor_factor=22/out_L)
                    compare_ax2 = feed2.plotSmithChart(ax=compare_ax2)
                    compare_ax2 = feed1.plotSmithChart(ax=compare_ax2)
            
            
                except Exception as e:
                    print(e)
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    print(exc_type, fname, exc_tb.tb_lineno)'''