#!/usr/bin/env python3
'''
This is code used to make the plots that were presented on 10/1/2020
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


if __name__ == '__main__':
    try:
        plt.close('all')
        cu_dir = '/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_copper_tab_testing_aug2020/s11/'
        cu_infiles = numpy.array(glob.glob(cu_dir + '*.csv'))
        cu_roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in cu_infiles])
        cu_linestyle = '-'
        
        al_dir = '/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_al_tube_testing_sep2020/s11/'
        al_infiles = numpy.array(glob.glob(al_dir + '*.csv'))
        al_roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in al_infiles])
        al_linestyle = '--'

        #Params
        linewidth = 3
        alpha = 0.9
        label_fontsize = 18
        legend_fontsize = 10
        plot_smith = True

        #PLOT 1, All data on same plot
        if plot_smith == True:
            all_fig = plt.figure()
            all_ax_lm = plt.subplot(2,1,1)
            all_ax_sm = get_smith(all_fig, rect = 212, plot_impedance = False, plot_ticks = True, plot_admittance = True, plot_labels = False)
            plt.ylabel('Im($\\Gamma$)',fontsize=label_fontsize)
            plt.xlabel('Re($\\Gamma$)',fontsize=label_fontsize)
            
        else:
            all_fig, all_ax_lm = plt.subplots(1,1)

        plt.tight_layout()
        all_ax_lm.set_xlabel('Freq (MHz)',fontsize=label_fontsize)
        all_ax_lm.set_ylabel('S11 (dB)',fontsize=label_fontsize)
        all_ax_lm.minorticks_on()
        all_ax_lm.grid(b=True, which='major', color='k', linestyle='-',alpha=0.75)
        all_ax_lm.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.25)
        
        all_roots = numpy.append(cu_roots, al_roots)
        for root_index, root in enumerate(all_roots):
            if root in cu_roots:
                linestyle = cu_linestyle
            else:
                linestyle = al_linestyle

            #Cap and inductor values
            shunt_inductor_value = root.split('/')[-1].split('_')[2] #nH
            if shunt_inductor_value.lower() == 'noshunt':
                shunt_inductor_value = 0
            else:
                shunt_inductor_value = float(shunt_inductor_value.lower().replace('nh',''))*1e-9 #H
            cap_value = float(root.split('/')[-1].split('_')[3].lower().replace('pf','').replace('p','.'))*1e-12 #F

            #Load data
            logmag_infile = root.replace('_FILLER','_LOGMAG')
            phase_infile = root.replace('_FILLER','_PHASE')
            
            freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
            freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)

            #Make feed, incase I want to plot Smith plot later.
            feed = SmithMatcher(freqs, logmag_vals, phase_vals, initial_capacitor_value=cap_value, initial_shunt_inductor_value=shunt_inductor_value, z_0=50)

            all_ax_lm.plot(freqs/1e6, logmag_vals,linestyle = linestyle,linewidth=linewidth,alpha=alpha,label='Series C\'s = %.1f pF, Shunt L = %i nH'%(cap_value*1e12,shunt_inductor_value*1e9))
            if plot_smith == True:
                all_ax_sm = feed.plotSmithChart(ax=all_ax_sm,linestyle=linestyle)

        all_ax_lm.set_xlim(0,1000)
        #Handle Legend
        #Get artists and labels for legend and chose which ones to display
        handles, labels = all_ax_lm.get_legend_handles_labels()
        display = numpy.arange(len(handles))#(0,len(cu_roots))

        #Create custom artists
        cuArtist = plt.Line2D((0,1),(0,0), color='k', linestyle=cu_linestyle,linewidth=linewidth)
        alArtist = plt.Line2D((0,1),(0,0), color='k', linestyle=al_linestyle,linewidth=linewidth)

        #Create legend from custom artist/label lists
        all_ax_lm.legend([handle for i,handle in enumerate(handles) if i in display]+[cuArtist,alArtist], [label for i,label in enumerate(labels) if i in display]+['Copper Tabbed Antenna', 'Aluminum Antenna'],fontsize=legend_fontsize,loc='lower left')





        #PLOT 2, Only Copper
        if plot_smith == True:
            cu_fig = plt.figure()
            cu_ax_lm = plt.subplot(2,1,1)
            cu_ax_sm = get_smith(cu_fig, rect = 212, plot_impedance = False, plot_ticks = True, plot_admittance = True, plot_labels = False)
            plt.ylabel('Im($\\Gamma$)',fontsize=label_fontsize)
            plt.xlabel('Re($\\Gamma$)',fontsize=label_fontsize)
            
        else:
            cu_fig, cu_ax_lm = plt.subplots(1,1)

        plt.tight_layout()
        cu_ax_lm.set_xlabel('Freq (MHz)',fontsize=label_fontsize)
        cu_ax_lm.set_ylabel('S11 (dB)',fontsize=label_fontsize)
        cu_ax_lm.minorticks_on()
        cu_ax_lm.grid(b=True, which='major', color='k', linestyle='-',alpha=0.75)
        cu_ax_lm.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.25)

        for root_index, root in enumerate(cu_roots):
            if root in cu_roots:
                linestyle = cu_linestyle
            else:
                linestyle = al_linestyle

            #Cap and inductor values
            shunt_inductor_value = root.split('/')[-1].split('_')[2] #nH
            if shunt_inductor_value.lower() == 'noshunt':
                shunt_inductor_value = 0
            else:
                shunt_inductor_value = float(shunt_inductor_value.lower().replace('nh',''))*1e-9 #H
            cap_value = float(root.split('/')[-1].split('_')[3].lower().replace('pf','').replace('p','.'))*1e-12 #F

            #Load data
            logmag_infile = root.replace('_FILLER','_LOGMAG')
            phase_infile = root.replace('_FILLER','_PHASE')
            
            freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
            freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)

            #Make feed, incase I want to plot Smith plot later.
            feed = SmithMatcher(freqs, logmag_vals, phase_vals, initial_capacitor_value=cap_value, initial_shunt_inductor_value=shunt_inductor_value, z_0=50)

            cu_ax_lm.plot(freqs/1e6, logmag_vals,linestyle = linestyle,linewidth=linewidth,alpha=alpha,label='Series C\'s = %.1f pF, Shunt L = %i nH'%(cap_value*1e12,shunt_inductor_value*1e9))
            if plot_smith == True:
                cu_ax_sm = feed.plotSmithChart(ax=cu_ax_sm,linestyle=linestyle)

        cu_ax_lm.set_xlim(0,1000)
        cu_ax_lm.legend(fontsize=legend_fontsize,loc='lower left')
        



        #PLOT 3, Only Al
        if plot_smith == True:
            al_fig = plt.figure()
            al_ax_lm = plt.subplot(2,1,1)
            al_ax_sm = get_smith(al_fig, rect = 212, plot_impedance = False, plot_ticks = True, plot_admittance = True, plot_labels = False)
            plt.ylabel('Im($\\Gamma$)',fontsize=label_fontsize)
            plt.xlabel('Re($\\Gamma$)',fontsize=label_fontsize)
            
        else:
            al_fig, al_ax_lm = plt.subplots(1,1)

        plt.tight_layout()
        al_ax_lm.set_xlabel('Freq (MHz)',fontsize=label_fontsize)
        al_ax_lm.set_ylabel('S11 (dB)',fontsize=label_fontsize)
        al_ax_lm.minorticks_on()
        al_ax_lm.grid(b=True, which='major', color='k', linestyle='-',alpha=0.75)
        al_ax_lm.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.25)

        for root_index, root in enumerate(al_roots):
            if root in cu_roots:
                linestyle = cu_linestyle
            else:
                linestyle = al_linestyle

            #Cap and inductor values
            antenna_name = root.split('/')[-1].split('_')[1]
            shunt_inductor_value = root.split('/')[-1].split('_')[2] #nH
            if shunt_inductor_value.lower() == 'noshunt':
                shunt_inductor_value = 0
            else:
                shunt_inductor_value = float(shunt_inductor_value.lower().replace('nh',''))*1e-9 #H
            cap_value = float(root.split('/')[-1].split('_')[3].lower().replace('pf','').replace('p','.'))*1e-12 #F

            # if shunt_inductor_value != 100*1e-9:
            #     continue

            #Load data
            logmag_infile = root.replace('_FILLER','_LOGMAG')
            phase_infile = root.replace('_FILLER','_PHASE')
            
            freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
            freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)

            #Make feed, incase I want to plot Smith plot later.
            feed = SmithMatcher(freqs, logmag_vals, phase_vals, initial_capacitor_value=cap_value, initial_shunt_inductor_value=shunt_inductor_value, z_0=50)

            al_ax_lm.plot(freqs/1e6, logmag_vals,linestyle = linestyle,linewidth=linewidth,alpha=alpha,label='Series C\'s = %.1f pF, Shunt L = %i nH, %s'%(cap_value*1e12,shunt_inductor_value*1e9,antenna_name))
            if plot_smith == True:
                al_ax_sm = feed.plotSmithChart(ax=al_ax_sm,linestyle=linestyle)

        al_ax_lm.set_xlim(0,1500)
        al_ax_lm.legend(fontsize=legend_fontsize,loc='lower left')





        #PLOT 4, Best Copper Only
        if plot_smith == True:
            best_cu_fig = plt.figure()
            best_cu_ax_lm = plt.subplot(2,1,1)
            best_cu_ax_sm = get_smith(best_cu_fig, rect = 212, plot_impedance = False, plot_ticks = True, plot_admittance = True, plot_labels = False)
            plt.ylabel('Im($\\Gamma$)',fontsize=label_fontsize)
            plt.xlabel('Re($\\Gamma$)',fontsize=label_fontsize)
            
        else:
            best_cu_fig, best_cu_ax_lm = plt.subplots(1,1)

        plt.tight_layout()
        best_cu_ax_lm.set_xlabel('Freq (MHz)',fontsize=label_fontsize)
        best_cu_ax_lm.set_ylabel('S11 (dB)',fontsize=label_fontsize)
        best_cu_ax_lm.minorticks_on()
        best_cu_ax_lm.grid(b=True, which='major', color='k', linestyle='-',alpha=0.75)
        best_cu_ax_lm.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.25)

        best_cu_roots = ['/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_copper_tab_testing_aug2020/s11/hpol_davebee_27nh_2p7pf_FILLER.csv',
                        '/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_copper_tab_testing_aug2020/s11/hpol_chuckape_27nh_2p7pf_FILLER.csv']
        for root_index, root in enumerate(best_cu_roots):
            if root in cu_roots:
                linestyle = cu_linestyle
            else:
                linestyle = al_linestyle

            #Cap and inductor values
            shunt_inductor_value = root.split('/')[-1].split('_')[2] #nH
            if shunt_inductor_value.lower() == 'noshunt':
                shunt_inductor_value = 0
            else:
                shunt_inductor_value = float(shunt_inductor_value.lower().replace('nh',''))*1e-9 #H
            cap_value = float(root.split('/')[-1].split('_')[3].lower().replace('pf','').replace('p','.'))*1e-12 #F

            #Load data
            logmag_infile = root.replace('_FILLER','_LOGMAG')
            phase_infile = root.replace('_FILLER','_PHASE')
            
            freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
            freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)

            #Make feed, incase I want to plot Smith plot later.
            feed = SmithMatcher(freqs, logmag_vals, phase_vals, initial_capacitor_value=cap_value, initial_shunt_inductor_value=shunt_inductor_value, z_0=50)

            best_cu_ax_lm.plot(freqs/1e6, logmag_vals,linestyle = linestyle,linewidth=linewidth,alpha=alpha,label='Series C\'s = %.1f pF, Shunt L = %i nH'%(cap_value*1e12,shunt_inductor_value*1e9))
            if plot_smith == True:
                best_cu_ax_sm = feed.plotSmithChart(ax=best_cu_ax_sm,linestyle=linestyle)

        best_cu_ax_lm.set_xlim(0,1000)
        best_cu_ax_lm.legend(fontsize=legend_fontsize,loc='lower left')
    except Exception as e:
        print('\nError in %s'%inspect.stack()[0][3])
        print(e)
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)