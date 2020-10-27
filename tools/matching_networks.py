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

from beacon.tools import field_fox as ff
from hpol_prototype.tools.pySmith import get_smith

import matplotlib.pyplot as plt
from matplotlib import lines


class SmithMatcher:
    '''
    Given initial logmag and unwrapped phase S11 measurements (with frequencies), this will provide the tools to warp
    the S11 based on Smith chart rules for shunt and in series components. This will use functions from the field fox
    class created in the BEACON analysis repository.  This assumes the specific part geometry of the hpol prototype.  
    The number of specified wings will determine the multiplicity of inductors added (1 shunt per element), as well as 
    the capacitor multiplicity (2 series per element).

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
        The impedance to be used for normalized impedance plots.  Typically set to 50 Ohms.s
    initial_capacitor_value : float
        The associated installed capacitors in the original logmag measurements.  This will be used to convert capacitor
        valus given to the actual part they need to be.  I.e. if you want to move the S11 a certain amount on the Smith
        chart, the difference in values of that component is what matters.  This should allow for inital input
        measurements to include shunt inductors, which might better extrapolate to different inductor values.
    initial_shunt_inductor_value : float
        The associated installed inductors in the original logmag measurements.  This will be used to convert inductor
        valus given to the actual part they need to be.  I.e. if you want to move the S11 a certain amount on the Smith
        chart, the difference in values of that component is what matters.  This should allow for inital input
        measurements to include shunt inductors, which might better extrapolate to different inductor values.
    readout_transformer_factor : float
            This fill be used to multiply z_0 when converting back to impedance from s11.  Can be used to simulate
            the effects of adding a transformer just before the sma. 
    '''
    def __init__(self, freqs, logmag, unwrapped_phase, initial_capacitor_value=2.7e-12, initial_shunt_inductor_value=0.0, z_0=50,readout_transformer_factor=1.0):
        try:
            self.n_elements = 1 #Theoretically we think this should be 3, but when trying to match data, 1 seems to work better.
            print('ASSUMING N_ELEMENTS = %i FOR THESE CALCULATIONS'%self.n_elements)
            self.n_series_caps_per_element = 2

            self.freqs = numpy.copy(freqs)
            self.logmag = numpy.copy(logmag)
            self.linmag = ff.logMagToLin(self.logmag)
            self.unwrapped_phase = unwrapped_phase
            re, im = ff.magPhaseToReIm(self.linmag,self.unwrapped_phase)
            self.complexs11 = re + 1.0j*im 
            self.z_0 = z_0
            self.readout_transformer_factor = readout_transformer_factor

            self.initial_capacitor_value = numpy.copy(initial_capacitor_value)
            self.initial_shunt_inductor_value = numpy.copy(initial_shunt_inductor_value)

            self.adjusted_capacitor_value = numpy.copy(initial_capacitor_value)
            self.adjusted_shunt_inductor_value = numpy.copy(initial_shunt_inductor_value)


            #These will be updated any time a component is added. 
            self.adjusted_logmag = numpy.copy(logmag)
            self.adjusted_linmag = ff.logMagToLin(self.adjusted_logmag)
            self.adjusted_unwrapped_phase = unwrapped_phase
            re, im = ff.magPhaseToReIm(self.adjusted_linmag,self.adjusted_unwrapped_phase)
            self.adjusted_complexs11 = re + 1.0j*im 

            self.calculateTransformedReadout()

            self.debug=False
        except Exception as e:
            print('\nError in %s'%inspect.stack()[0][3])
            print(e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)

    def reset(self):
        '''
        This will set all adjusted values back to the original initiated values.
        '''
        self.adjusted_logmag = numpy.copy(self.logmag)
        self.adjusted_linmag = numpy.copy(self.linmag)
        self.adjusted_unwrapped_phase = numpy.copy(self.unwrapped_phase)
        self.adjusted_complexs11 = numpy.copy(self.complexs11)

        self.adjusted_capacitor_value = numpy.copy(self.initial_capacitor_value)
        self.adjusted_shunt_inductor_value = numpy.copy(self.initial_shunt_inductor_value)

    def calculateTransformedReadout(self):
        try:
            #Read in impedance as it IS (with z_0 = 50 Ohm)
            z = ff.logToComplexZ(self.adjusted_logmag, self.adjusted_unwrapped_phase , z_0=self.z_0)

            #Read out impedance as seen by SMA after transformer (if sma has transformer multiplying it's imepdance)
            self.transformed_adjusted_complexs11 = (z - self.readout_transformer_factor*self.z_0)/(z + self.readout_transformer_factor*self.z_0)
            self.transformed_adjusted_linmag = numpy.abs(self.transformed_adjusted_complexs11)
            self.transformed_adjusted_logmag = ff.linToLogMag(self.transformed_adjusted_linmag)
            self.transformed_adjusted_unwrapped_phase = numpy.unwrap(numpy.angle(self.transformed_adjusted_complexs11,deg=True),discont=360.0)
        except Exception as e:
            print('\nError in %s'%inspect.stack()[0][3])
            print(e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)

    def forceAdjustedAsInitial(self):
        '''
        This will set the current adjusted as the "original" input response and cap and inductor values.  T
        his change cannot be undone, as it overwrites the normal touchpoint for self.reset.  This could be helpful if 
        you wanted to test the effects of varying values with an assumed additional offset (for instance artificially 
        adding resistance), but should be used with caution.
        '''
        print('Warning:  The initial S11 for this feed has been overwritten with current adjusted values.')
        self.logmag = numpy.copy(self.adjusted_logmag)
        self.linmag = numpy.copy(self.adjusted_linmag)
        self.unwrapped_phase = numpy.copy(self.adjusted_unwrapped_phase)
        self.complexs11 = numpy.copy(self.adjusted_complexs11)

        self.initial_capacitor_value = numpy.copy(self.adjusted_capacitor_value)
        self.initial_shunt_inductor_value = numpy.copy(self.adjusted_shunt_inductor_value)

    def getPercentBelow(self, db=-5,freq_low_cut_MHz=100,freq_high_cut_MHz=800):
        '''
        This will return the percentage of points in the band that are below the given db value.  This will be 
        calculated using the adjusted S11 logmag curve.
        '''
        cut = numpy.logical_and(self.freqs >= freq_low_cut_MHz, self.freqs <= freq_high_cut_MHz)
        return sum(self.adjusted_logmag < db)/len(self.adjusted_logmag)

    def plotCurrentLogMagS11(self,ax=None,fontsize=16,leg_fontsize=14,label='',plot_cut_ll=0,plot_cut_ul=1500,plot_transformed=False):
        try:
            plot_cut = numpy.logical_and(self.freqs/1e6 > plot_cut_ll, self.freqs/1e6 <= plot_cut_ul)

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

            if plot_transformed == True:
                _ax.plot(self.freqs[plot_cut]/1e6, self.transformed_adjusted_logmag[plot_cut],label=label,alpha=1.0,linestyle='-')#,color=color)
            else:
                _ax.plot(self.freqs[plot_cut]/1e6, self.adjusted_logmag[plot_cut],label=label,alpha=1.0,linestyle='-')#,color=color)
            _ax.legend(loc='lower right')

            return _ax
        except Exception as e:
            print('\nError in %s'%inspect.stack()[0][3])
            print(e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)

    def plotSmithChart(self,ax=None,linestyle='-',fontsize=16,leg_fontsize=14,label='',plot_cut_ll=0,plot_cut_ul=1500,plot_transformed=False):
        try:
            plot_cut = numpy.logical_and(self.freqs/1e6 > plot_cut_ll, self.freqs/1e6 <= plot_cut_ul)
            
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

            if plot_transformed == True:
                _ax.plot(numpy.real(self.transformed_adjusted_complexs11)[plot_cut], numpy.imag(self.transformed_adjusted_complexs11)[plot_cut],label=label,alpha=1.0,linestyle=linestyle)#,color=color)
            else:
                _ax.plot(numpy.real(self.adjusted_complexs11)[plot_cut], numpy.imag(self.adjusted_complexs11)[plot_cut],label=label,alpha=1.0,linestyle=linestyle)#,color=color)
            _ax.legend(loc='upper right')

            return _ax
        except Exception as e:
            print('\nError in %s'%inspect.stack()[0][3])
            print(e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)

    def swapRLC(self, part, val, verbose=False):
        '''
        This will shift the reflection coefficient from the existing most recent adjusted logmag, phase, to the new
        appropriate value.  This accounts for SWAPPING a part, which means it accounts for the existing rlc values
        for the calculation of shifting along the curve.  For parts that are not currently
        on the feed (pc, sr, pr, sl), this will not use any assumed existing part, and should only be used for tests.  
        These components shouldn't be used without caution. They will also not be added on a per-element basis as their
        multiplicity is unknown. 
    
        Part should be in 'sc, pc, sr, pr, sl, pl'
        '''
        #Z = R + iX
        #Y = 1/Z = G + iB
        try:
            z = ff.logToComplexZ(self.adjusted_logmag, self.adjusted_unwrapped_phase )
            norm_z = z/self.z_0
            if part == 'sr':
                #R += val/self.z_0
                #Series resistor
                if verbose:
                    print('WARNING!!! Adding Series resistor')
                norm_z += val/self.z_0
                out_z = norm_z*self.z_0
            elif part == 'pr':
                #Parallel resistor
                if verbose:
                    print('WARNING!!! Adding Parallel resistor')
                norm_y = 1.0/norm_z
                norm_y += self.z_0/val
                norm_z = 1.0/norm_y
                out_z = norm_z*self.z_0
            elif part == 'sc':
                #Series capacitor
                old_sc = self.n_elements*self.adjusted_capacitor_value/self.n_series_caps_per_element
                new_sc = self.n_elements*val/self.n_series_caps_per_element
                if verbose:
                    print('Adding Series capacitor')
                norm_z -= 1.0j * (1.0/new_sc - 1.0/old_sc) /(2.0*numpy.pi*self.freqs*self.z_0)
                out_z = norm_z*self.z_0
                self.adjusted_capacitor_value = val
            elif part == 'pc':
                #Parallel capacitor
                if verbose:
                    print('WARNING!!! Adding Parallel capacitor')
                norm_y = 1.0/norm_z
                norm_y += 1.0j * (2.0*numpy.pi*self.freqs*val*self.z_0)
                norm_z = 1.0/norm_y
                out_z = norm_z*self.z_0
            elif part == 'sl':
                #Series inductor
                if verbose:
                    print('WARNING!!! Adding Series inductor')
                norm_z += 1.0j * (2.0*numpy.pi*self.freqs*val)/(self.z_0)
                out_z = norm_z*self.z_0
            elif part == 'pl':
                #Parallel inductor
                old_pl = self.adjusted_shunt_inductor_value/self.n_elements
                new_pl = val/self.n_elements
                if verbose:
                    print('Adding Parallel inductor')
                norm_y = 1.0/norm_z
                if val == 0:
                    print('Assuming no shunt component')
                    norm_y += 1.0j * (1.0/old_pl) * self.z_0/(2.0*numpy.pi*self.freqs) #Should allow to predict no shunt case
                elif old_pl == 0.0:
                    norm_y -= 1.0j * (1.0/new_pl) * self.z_0/(2.0*numpy.pi*self.freqs) #What I was doing before
                else:
                    norm_y -= 1.0j * (1.0/new_pl - 1.0/old_pl) * self.z_0/(2.0*numpy.pi*self.freqs) #Should handle everything else.
                norm_z = 1.0/norm_y
                out_z = norm_z*self.z_0
                self.adjusted_shunt_inductor_value = val
            else:
                print('COMPONENT NOT ADDED, VALUE NOT ACCEPTED.')

            self.adjusted_complexs11 = (out_z - self.z_0)/(out_z + self.z_0)
            self.adjusted_linmag = numpy.abs(self.adjusted_complexs11)
            self.adjusted_logmag = ff.linToLogMag(self.adjusted_linmag)
            self.adjusted_unwrapped_phase = numpy.unwrap(numpy.angle(self.adjusted_complexs11,deg=True),discont=360.0)

            self.calculateTransformedReadout()

        except Exception as e:
            print('\nError in %s'%inspect.stack()[0][3])
            print(e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)

    def addRLC(self, part, val, verbose=False):
        '''
        This will shift the reflection coefficient from the existing most recent adjusted logmag, phase, to the new
        appropriate value.  This does not account for SWAPPING a part, and will simply blindly add to the complex
        impedance without considering existing components.  

        Part should be in 'sc, pc, sr, pr, sl, pl'
        '''
        #Z = R + iX
        #Y = 1/Z = G + iB
        try:
            z = ff.logToComplexZ(self.adjusted_logmag, self.adjusted_unwrapped_phase )
            norm_z = z/self.z_0
            if part == 'sr':
                #R += val/self.z_0
                #Series resistor
                if verbose:
                    print('Adding Series resistor')
                norm_z += val/self.z_0
                out_z = norm_z*self.z_0
            elif part == 'pr':
                #Parallel resistor
                if verbose:
                    print('Adding Parallel resistor')
                norm_y = 1.0/norm_z
                norm_y += self.z_0/val
                norm_z = 1.0/norm_y
                out_z = norm_z*self.z_0
            elif part == 'sc':
                #Series capacitor
                if verbose:
                    print('Adding Series capacitor')
                norm_z -= 1.0j /(2.0*numpy.pi*self.freqs*val*self.z_0)
                out_z = norm_z*self.z_0
            elif part == 'pc':
                #Parallel capacitor
                if verbose:
                    print('Adding Parallel capacitor')
                norm_y = 1.0/norm_z
                norm_y += 1.0j * (2.0*numpy.pi*self.freqs*val*self.z_0)
                norm_z = 1.0/norm_y
                out_z = norm_z*self.z_0
            elif part == 'sl':
                #Series inductor
                if verbose:
                    print('Adding Series inductor')
                norm_z += 1.0j * (2.0*numpy.pi*self.freqs*val)/(self.z_0)
                out_z = norm_z*self.z_0
            elif part == 'pl':
                #Parallel inductor
                if verbose:
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
            self.adjusted_unwrapped_phase = numpy.unwrap(numpy.angle(self.adjusted_complexs11,deg=True),discont=360.0)

            self.calculateTransformedReadout()
        except Exception as e:
            print('\nError in %s'%inspect.stack()[0][3])
            print(e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)



class HpolTriWingFeed(SmithMatcher):
    '''
    THIS IS DEPRECATED, DO NOT USE.  SMITH MATCHER SHOULD BE USED INSTEAD!!!

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
    initial_capacitor_value : float
        The associated installed capacitors in the original logmag measurements.  This will be used to convert capacitor
        valus given to the actual part they need to be.  I.e. if you want to move the S11 a certain amount on the Smith
        chart, the difference in values of that component is what matters.  This should allow for inital input
        measurements to include shunt inductors, which might better extrapolate to different inductor values.
    initial_shunt_inductor_value : float
        The associated installed inductors in the original logmag measurements.  This will be used to convert inductor
        valus given to the actual part they need to be.  I.e. if you want to move the S11 a certain amount on the Smith
        chart, the difference in values of that component is what matters.  This should allow for inital input
        measurements to include shunt inductors, which might better extrapolate to different inductor values.
    '''
    def __init__(self, freqs, logmag, unwrapped_phase, initial_capacitor_value=2.7e-12, initial_shunt_inductor_value=0.0, z_0=50):
        try:
            print('WARNING!!! HPOLTRIWINGFEED is DEPRECATED, DO NOT USE.')
            super().__init__(freqs, logmag, unwrapped_phase, z_0=z_0)
            self.initial_capacitor_value = initial_capacitor_value
            self.initial_shunt_inductor_value = initial_shunt_inductor_value
            self.n_elements = 3
            self.n_series_caps_per_element = 2

        except Exception as e:
            print('\nError in %s'%inspect.stack()[0][3])
            print(e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)

    def convertAddedComponentToSwappedPart(self, part, val, verbose=False):
        '''
        addRLC will simply adjust the Smith chart as expected for adding that quantity/part.  In the lab we often
        don't just add new parts on top, we tweak the value of the existing part by swapping it out.  This does the
        conversion from the value of a part used to shift the Smith chart in addRLC, to the part that should be 
        placed.  This assumes the effect circuit geometry/multiplicity of parts of the HpolTriWingFeed.
        '''
        print('I want to implement this and try to make the dB saturation plots.')


    def addShuntInductor(self, val, resistance=0.0):
        '''
        This will do the math to emulate adding 3 shunt inductors of the given value (assumed usits of H), one inductor
        per wing of the center feed.

        The kwarg for resistance allows you to attempt to account for the added resistance from the inductor as 
        described on the spec sheet.  Beacuse we are talking about a resistance intrinsic to the inductor, it is also
        assumed to be added in parallel to the load.  
        ''' 
        for wing in range(self.n_elements):
            self.addRLC('pl',val)
            '''
            if resistance != 0:
                print('Still unsure if this is the right call, seems to short the antenna.')
                self.addRLC('sr',resistance)
            '''
    def addSeriesCapacitor(self, val):
        '''
        This will do the math to emulate adding adding additional capacitance to the network, assuming you are adding 6
        capacitors, in 3 parallel groups of 2 in series capacitors.
        ''' 
        for wing in range(self.n_elements):
            self.addRLC('sc',val/self.n_series_caps_per_element)

    def getPercentBelow(self, db=-5):
        '''
        This will return the percentage of points in the band that are below the given db value.  This will be 
        calculated using the adjusted S11 logmag curve.
        '''
        return sum(self.adjusted_logmag < db)/len(self.adjusted_logmag)

    def add22nHShuntInductor(self,effective_inductor_factor=1.0):
        #Z = R + iX
        #Y = 1/Z = G + iB
        #pulling data from pg 2 of https://www.mouser.com/datasheet/2/54/ce201210-777343.pdf
        try:
            freqs = 1e6*numpy.array([10.14251747516137, 10.831709055888476, 11.567731715399361, 12.353767660211613, 13.193215330134526, 14.089704091481606, 15.047109928697665, 16.06957220224368, 17.161511545192308, 18.327648975910368, 19.573026309462186, 20.903027955983063, 22.32340420026828, 23.840296063227406, 25.460261852692966, 27.19030551837611, 29.037906933562674, 31.011054234472823, 33.1182783571042, 35.368689920879184, 37.77201861856329, 40.338655282757664, 43.07969681084113, 46.00699414259512, 49.133203497943896, 52.471841096337165, 56.0373415943554, 59.845120494193495, 63.91164079284757, 68.25448416016302, 72.89242695348393, 77.84552139755414, 83.13518228065288, 88.78427954179689, 94.81723714931186, 101.26013869827587, 108.14084018338724, 115.48909043483285, 123.3366597378625, 131.71747719215892, 140.66777740487868, 150.2262571515941, 160.43424268246218, 171.3358683969718, 182.97826765977254, 195.41177658257968, 208.69015165321247, 222.87080215268514, 238.01503836521218, 254.18833665426774, 271.4606225507593, 289.9065730772528, 309.60593961535295, 330.6438927121616, 353.1113903165897, 377.10557103760095, 402.73017412464344, 430.0959879860635, 459.321329184678, 490.53255398145194, 523.8646046389445, 559.4615928464851, 597.4774227895181, 638.0764565569717, 681.4342247635534, 727.7381854593541, 777.1885346079248, 829.9990716369375, 886.3981238036572, 946.6295333717304, 1007.9372333072872, 1060.4636727717857, 1107.7542495737584, 1125.774615744299, 1170.3685102199822, 1199.1963227759618, 1231.359699623225, 1268.7111706000474, 1289.735247936239, 1299.4064736146222, 1315.0315039760671, 1334.3841452242648, 1349.1492920415396, 1370.3033524781254, 1387.7021067008495, 1396.0205559996566])
            #inductor_values = 1e-9*numpy.array([25.92922758914868, 25.88704077719219, 26.003220006504996, 25.93978502990999, 25.88704077719219, 25.88704077719219, 25.82388922036946, 25.63535739861822, 25.54160827109554, 25.35513729080619, 24.986269483408535, 24.91517080873166, 24.864509939178284, 24.864509939178284, 24.854390133233892, 24.662894473345844, 24.64282306098564, 24.80385285110829, 25.047372665954132, 25.088191112904397, 25.006620630552003, 24.976100121572536, 24.70308635546962, 24.64282306098564, 24.622767983350247, 24.53272192002691, 24.53272192002691, 24.433056903015306, 24.423112697595954, 24.294205729953607, 23.999318858437093, 24.038429344663733, 24.20536122680803, 24.225076329469314, 24.107026123744223, 23.746647126944268, 23.804718849704127, 24.009090510380375, 24.264554721323336, 24.313993195860526, 24.313993195860526, 24.313993195860526, 24.225076329469314, 23.911552763776044, 23.679075809830834, 23.698362256587014, 24.009090510380375, 24.20536122680803, 24.313993195860526, 24.313993195860526, 24.313993195860526, 24.383376332167277, 24.37345234664115, 24.18566216889591, 24.37345234664115, 24.723206851874743, 24.965934898646516, 25.12907607963034, 25.30358182234648, 25.76089172206438, 26.365626247892514, 26.82028546298194, 27.606776333048668, 29.142884078430214, 29.972838669250713, 30.57564151198693, 31.999340784935402, 33.64628186540018, 36.05719040778525, 39.403314059702524, 42.33161261194575, 47.1408516508453, 51.98472167419724, 57.238370291719946, 63.051187283834835, 70.53788959903306, 78.52079325372053, 90.49662568681987, 102.58517888493259, 117.07224270242266, 134.30497459460125, 153.1571571481011, 178.3437385455982, 203.68136870353183, 229.6126372841107, 254.5214610474032])
            inductor_values = numpy.ones_like(freqs)*22e-9
            inductor_curve = scipy.interpolate.interp1d(freqs, inductor_values, kind='cubic', assume_sorted=False,fill_value='extrapolate')

            interpolated_inductor_values = (1/effective_inductor_factor)*inductor_curve(self.freqs)

            z = ff.logToComplexZ(self.adjusted_logmag, self.adjusted_unwrapped_phase )
            norm_z = z/self.z_0
            #Parallel inductor

            norm_y = 1.0/norm_z
            norm_y -= 1.0j * self.z_0/(2.0*numpy.pi*self.freqs*interpolated_inductor_values)
            norm_z = 1.0/norm_y
            out_z = norm_z*self.z_0


            self.adjusted_complexs11 = (out_z - self.z_0)/(out_z + self.z_0)

            self.adjusted_linmag = numpy.abs(self.adjusted_complexs11)
            self.adjusted_logmag = ff.linToLogMag(self.adjusted_linmag)
            self.adjusted_unwrapped_phase = numpy.unwrap(numpy.angle(self.adjusted_complexs11,deg=True),discont=360.0)

        except Exception as e:
            print('\nError in %s'%inspect.stack()[0][3])
            print(e)
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)


if __name__ == '__main__':
    plt.close('all')
    datapath ='/home/dsouthall/Projects/Greenland/hpol_prototype/data/hpol_copper_tab_testing_aug2020/s11/'
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

    #sm = HpolTriWingFeed(freqs, logmag_vals, phase_vals, z_0=50)
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