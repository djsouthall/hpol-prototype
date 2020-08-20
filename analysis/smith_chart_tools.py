#!/usr/bin/env python3
'''
This is meant to contain a class that, once give an S11 curve in logmag and phase, can produce a Smith chart for that
curve, and also predict how that smith chart would vary with different components. 
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
from pySmith import get_smith

import matplotlib.pyplot as plt
from matplotlib import lines


class SmithMatcher:
    '''
    Given initial logmag and unwrapped phase S11 measurements (with frequencies), this will provide the tools to warp
    the S11 based on Smith chart rules for shunt and in series components. This will use functions from the field fox
    class created in the BEACON analysis repository.  

    Parameters
    ----------
    freqs : numpy.array of floats
        The frequencies as read out by the field fox.  Should be given in Hz.
    logmag : numpy.array of floats
        The log magnitudes of the S11 measurement as readout by the field fox.  Should be given in dB as measured in a
        field fox network analyzer.
    unwrapped_phase : numpy.array of floats
        The unwrapped phase for the S11 measurement. Should be given in degrees as measured in a field fox network
        analyzer.
    z_0 : float
        The impedance to be used for normalized impedance plots.  Typically set to 50 Ohms.

    '''
    def __init__(self, freqs, logmag, unwrapped_phase, z_0=50):
        try:
            self.freqs = freqs
            self.logmag = logmag
            self.linmag = ff.logMagToLin(self.logmag)
            self.unwrapped_phase = unwrapped_phase
            re, im = ff.magPhaseToReIm(self.linmag,self.unwrapped_phase)
            self.complexs11 = re + 1.0j*im 
            self.z_0 = z_0


            #These will be updated any time a component is added. 
            self.adjusted_logmag = logmag
            self.adjusted_linmag = ff.logMagToLin(self.adjusted_logmag)
            self.adjusted_unwrapped_phase = unwrapped_phase
            self.unwrapped_phase = unwrapped_phase
            re, im = ff.magPhaseToReIm(self.adjusted_linmag,self.adjusted_unwrapped_phase)
            self.adjusted_complexs11 = re + 1.0j*im 


        except Exception as e:
            print('\nError in %s'%inspect.stack()[0][3])
            print(e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)

    def plotCurrentLogMagS11(self,ax=None,fontsize=16,leg_fontsize=14,label=''):
        try:
            plot_cut_ll = 100            
            plot_cut = self.freqs/1e6 > plot_cut_ll

            if ax is None:
                #PLOT PREPPING
                _fig = plt.figure()
                _ax = plt.subplot(1,1,1)
                plt.ylabel('dB',fontsize=fontsize)
                plt.xlabel('MHz',fontsize=fontsize)
                plt.minorticks_on()
                plt.grid(b=True, which='major', color='k', linestyle='-')
                plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)
            else:
                #Use passed ax
                _ax = ax

            _ax.plot(self.freqs[plot_cut]/1e6, self.adjusted_logmag[plot_cut],label=label,alpha=1.0,linestyle='-')#,color=color)
            _ax.legend(loc='lower right')

            return _ax
        except Exception as e:
            print('\nError in %s'%inspect.stack()[0][3])
            print(e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)

    def plotSmithChart(self,ax=None,fontsize=16,leg_fontsize=14,label=''):
        try:
            if ax is None:
                #PLOT PREPPING
                _fig = plt.figure()
                _ax = get_smith(_fig, rect = 111, plot_impedance = False, plot_ticks = True, plot_admittance = True, plot_labels = False)
                plt.ylabel('Im($\\Gamma$)',fontsize=fontsize)
                plt.xlabel('Re($\\Gamma$)',fontsize=fontsize)
                plt.minorticks_on()
                plt.grid(b=True, which='major', color='k', linestyle='-')
                plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)
            else:
                #Use passed ax
                _ax = ax
                print('Using predefined ax')

            _ax.plot(numpy.real(self.adjusted_complexs11)[plot_cut], numpy.imag(self.adjusted_complexs11)[plot_cut],label=label,alpha=1.0,linestyle='-')#,color=color)
            _ax.legend(loc='upper right')

            return _ax
        except Exception as e:
            print('\nError in %s'%inspect.stack()[0][3])
            print(e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)

    def addRLC(self, part, val):
        '''
        This will shift the reflection coefficient from the existing most recent adjusted logmag, phase, to the new
        appropriate value. 
    
        part should be in 'sc, pc, sr, pr, sl, pl'
        '''
        #Z = R + iX
        #Y = 1/Z = G + iB
        z = ff.logToComplexZ(self.adjusted_logmag, self.unwrapped_phase, z_0 = self.z_0 )
        norm_z = z/self.z_0
        if part == 'sr':
            #R += val/self.z_0
            #Series resistor
            print('Adding Series resistor')
            norm_z += val/self.z_0
            out_z = norm_z*self.z_0
        elif part == 'pr':
            #Parallel resistor
            print('Adding Parallel resistor')
            norm_y = 1.0/norm_z
            norm_y += self.z_0/val
            norm_z = 1.0/norm_y
            out_z = norm_z*self.z_0
        elif part == 'sc':
            #Series capacitor
            print('Adding Series capacitor')
            norm_z -= 1.0j /(2.0*numpy.pi*self.freqs*val*self.z_0)
            out_z = norm_z*self.z_0
        elif part == 'pc':
            #Parallel capacitor
            print('Adding Parallel capacitor')
            norm_y = 1.0/norm_z
            norm_y += 1.0j * (2.0*numpy.pi*self.freqs*val*self.z_0)
            norm_z = 1.0/norm_y
            out_z = norm_z*self.z_0
        elif part == 'sl':
            #Series inductor
            print('Adding Series inductor')
            norm_z += 1.0j * (2.0*numpy.pi*self.freqs*val)/(self.z_0)
            out_z = norm_z*self.z_0
        elif part == 'pl':
            #Parallel inductor
            print('Adding Parallel inductor')
            norm_y = 1.0/norm_z
            norm_y -= 1.0j * self.z_0/(2.0*numpy.pi*self.freqs*val)
            norm_z = 1.0/norm_y
            out_z = norm_z*self.z_0
        else:
            print('COMPONENT NOT ADDED, VALUE NOT ACCEPTED.')

        self.adjusted_complexs11 = (out_z - self.z_0)/(out_z + self.z_0)

        self.adjusted_linmag = numpy.abs(self.adjusted_complexs11)
        self.adjusted_logmag = ff.linToLogMag(self.adjusted_linmag)
        self.unwrapped_phase = numpy.unwrap(numpy.angle(self.adjusted_complexs11))




class HpolTriWingFeed(SmithMatcher):
    '''
    This is a specific implementation of the Smith Matcher that adds function specific to the Tri-Wing centerfeed we are
    testing.  Specifically it adds functions to convert the realworld shunt inductor and series capacitor values to
    their equivalent values for changing Smith Charts.  It may also expand scope as I see fit.

    Parameters
    ----------
    freqs : numpy.array of floats
        The frequencies as read out by the field fox.  Should be given in Hz.
    logmag : numpy.array of floats
        The log magnitudes of the S11 measurement as readout by the field fox.  Should be given in dB as measured in a
        field fox network analyzer.
    unwrapped_phase : numpy.array of floats
        The unwrapped phase for the S11 measurement. Should be given in degrees as measured in a field fox network
        analyzer.
    z_0 : float
        The impedance to be used for normalized impedance plots.  Typically set to 50 Ohms.

    '''
    def __init__(self, freqs, logmag, unwrapped_phase, z_0=50):
        try:
            super().__init__(freqs, logmag, unwrapped_phase, z_0=z_0)
        except Exception as e:
            print('\nError in %s'%inspect.stack()[0][3])
            print(e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)

    def addShuntInductor(self, val):
        '''
        This will do the math to emulate adding 3 shunt inductors of the given value (assumed usits of H), one inductor
        per wing of the center feed.
        ''' 
        self.addRLC('pl',val/3.0)
    def addSeriesCapacitor(self, val):
        '''
        This will do the math to emulate adding adding additional capacitance to the network, assuming you are adding 6
        capacitors, in 3 parallel groups of 2 in series capacitors.
        ''' 
        self.addRLC('pl',3.0*val/2.0)


if __name__ == '__main__':
    plt.close('all')
    datapath ='/home/dsouthall/Projects/Greenland/hpol-prototype/data/hpol_copper_tab_testing_aug2020/s11/'
    infiles = numpy.array(glob.glob(datapath + '*.csv'))
    roots = numpy.unique([infile.replace('_PHASE','_FILLER').replace('_LOGMAG','_FILLER') for infile in infiles])

    #PLOT PREPPING
    fontsize=16
    leg_fontsize=14
    #PLOT 1, LOGMAG
    logmag_fig = plt.figure()
    logmag_ax = plt.subplot(1,1,1)
    plt.ylabel('dB',fontsize=fontsize)
    plt.xlabel('MHz',fontsize=fontsize)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='k', linestyle='-')
    plt.grid(b=True, which='minor', color='tab:gray', linestyle='--',alpha=0.5)

    for root in roots:
        logmag_infile = root.replace('_FILLER','_LOGMAG')
        phase_infile = root.replace('_FILLER','_PHASE')

        label = logmag_infile.split('/')[-1].replace('.csv','').replace('_PHASE','').replace('_LOGMAG','').title().replace('pf',' pF').title().replace('Pf',' pF').replace('Nh',' nH').replace('nh',' nH').replace('2p7','2.7').replace('_', ' ').replace('NOSHUNT','No Shunt')

        freqs, logmag_vals = ff.readerFieldFox(logmag_infile,header=17)
        freqs, phase_vals = ff.readerFieldFox(phase_infile,header=17)
        
        plot_cut_ll = 100            
        plot_cut = freqs/1e6 > plot_cut_ll

        logmag_ax.plot(freqs[plot_cut]/1e6, logmag_vals[plot_cut],label=label,alpha=1.0,linestyle='-')#,color=color)

    plt.legend(loc='lower right')


    #Testing class

    sm = HpolTriWingFeed(freqs, logmag_vals, phase_vals, z_0=50)
    sm.plotCurrentLogMagS11(label='test')
    part_test_dict = {'sr':10,'pr':10,'sc':2.7e-12,'pc':2.7e-12,'sl':100e-9,'pl':100e-9}

    for key in list(part_test_dict.keys()):
        sm = HpolTriWingFeed(freqs, logmag_vals, phase_vals, z_0=50)
        ax_1 = sm.plotCurrentLogMagS11(ax=None,fontsize=16,leg_fontsize=14,label=label)
        ax_2 = sm.plotSmithChart(ax=None,fontsize=16,leg_fontsize=14,label=label)
        original = sm.adjusted_complexs11
        sm.addRLC(key, part_test_dict[key])
        ax_1 = sm.plotCurrentLogMagS11(ax=ax_1,fontsize=16,leg_fontsize=14,label='Added %s with val %0.3e'%(key, part_test_dict[key]))
        ax_2 = sm.plotSmithChart(ax=ax_2,fontsize=16,leg_fontsize=14,label='Added %s with val %0.3e'%(key, part_test_dict[key]))
        adjusted = sm.adjusted_complexs11
        print(numpy.all(original == adjusted))


    plt.figure()
    plt.plot(sm.freqs, sm.unwrapped_phase)
    plt.plot(sm.freqs, sm.adjusted_unwrapped_phase)

    sm = HpolTriWingFeed(freqs, logmag_vals, phase_vals, z_0=50)
    ax_1 = sm.plotCurrentLogMagS11(ax=None,fontsize=16,leg_fontsize=14,label=label)
    ax_2 = sm.plotSmithChart(ax=None,fontsize=16,leg_fontsize=14,label=label)
    for inductor_value in [47e-9,56e-9,100e-9,110e-9]:
        sm = HpolTriWingFeed(freqs, logmag_vals, phase_vals, z_0=50) #resetting value, should add function for this.
        sm.addShuntInductor(inductor_value)
        ax_1 = sm.plotCurrentLogMagS11(ax=ax_1,fontsize=16,leg_fontsize=14,label='Added Shunt Inductors with value = %0.2f nH'%(inductor_value*1e9))
        ax_2 = sm.plotSmithChart(ax=ax_2,fontsize=16,leg_fontsize=14,label='Added Shunt Inductors with value = %0.2f nH'%(inductor_value*1e9))
        sm.addShuntInductor(-inductor_value)