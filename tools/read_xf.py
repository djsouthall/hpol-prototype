# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 15:50:59 2020

@author: Originally by Bryan, Edited by Dan
"""

import pandas as pd
import numpy  as np
import matplotlib.pyplot as plt
plt.ion()


file_S11 = r'/home/dsouthall/Projects/Greenland/hpol_prototype/data/xf/sample/600mm 8in OD Closed Al S11.csv'
file_Z = r'/home/dsouthall/Projects/Greenland/hpol_prototype/data/xf/sample/600mm 8in OD Closed Al Impedance.csv'

#cutoff data beyond range of frequency interest given by freq_cutoff
freq_cutoff = 1.0
axes_label_fontsize = 15.0
tick_fontsize = 14.0
legend_fontsize = 20.0
title_fontsize = 25.0

fig = plt.figure(figsize=(14.0,10.0))
ax = fig.add_subplot(1, 1, 1)
plt.minorticks_on()
major_x_ticks = np.arange(0, 1.1, 0.2)
minor_x_ticks = np.arange(0, 1.1, 0.05)
major_y_ticks = np.arange(0, -100.0, -5.0)
minor_y_ticks = np.arange(0, -100.0, -1.0)

ax.set_xticks(major_x_ticks)
ax.set_xticks(minor_x_ticks, minor=True)
ax.set_yticks(major_y_ticks)
ax.set_yticks(minor_y_ticks, minor=True)
ax.grid(True, which = 'major', linestyle='-', color='black')
ax.grid(True, which = 'minor', linestyle='--')
#min_S11 = -10

xf_S11_df = pd.read_csv(file_S11)
xf_S11_df.drop(xf_S11_df.index[0], inplace = True)
xf_S11_df.rename(columns={"FFT Dimension (GHz)": "GHz"}, inplace= True)
xf_S11_df = xf_S11_df[xf_S11_df['GHz']<=freq_cutoff]
xf_S11_df['Re squared'] = xf_S11_df['Re( S11 )'].astype(float).apply(lambda x: x**2)
xf_S11_df['Im squared'] = xf_S11_df['Im( S11 )'].astype(float).apply(lambda x: x**2)
xf_S11_df['S11 (mag)'] = xf_S11_df['Re squared']+xf_S11_df['Im squared']
xf_S11_df['S11 (mag)'] = xf_S11_df['S11 (mag)'].astype(float).apply(lambda x: 10*np.log10(x))


plt.xticks(fontsize = tick_fontsize)
plt.yticks(fontsize = tick_fontsize)
plt.xlabel('Freq (GHz)', fontsize = axes_label_fontsize)
plt.ylabel('S11 (dB)', fontsize = axes_label_fontsize)
plt.title("S11 (dB) vs Freq (GHz)", fontsize = title_fontsize)
#plt.legend(loc="center left", fontsize = legend_fontsize)
plt.plot(xf_S11_df['GHz'], xf_S11_df['S11 (mag)'])
#plt.savefig("S11.csv")

# =============================================================================

#import csv and skip first row since XFdtd first row is broken
xf_Z_df = pd.read_csv(file_Z, skiprows = [1])
xf_Z_df.rename(columns={"FFT Dimension (GHz)": "GHz","Re( Z ) (ohm)": "Re(Z)", "Im( Z ) (ohm)": "Im(Z)"}, inplace=True)
xf_Z_df = xf_Z_df[xf_Z_df['GHz']<=freq_cutoff]

#create numpy arrays from dataframe
freq_GHz = xf_Z_df['GHz']
Z_Re_antenna = xf_Z_df['Re(Z)']
Z_Im_antenna = xf_Z_df['Im(Z)']

fig = plt.figure(figsize=(14.0,10.0))
ax = fig.add_subplot(1, 1, 1)
plt.minorticks_on()
major_x_ticks = np.arange(0, freq_cutoff + 0.1, 0.2)
minor_x_ticks = np.arange(0, freq_cutoff + 0.1, 0.05)
major_y_ticks = np.arange(-1600, 1600.0, 200.0)
minor_y_ticks = np.arange(-1600, 1600.0, 50.0)

ax.set_xticks(major_x_ticks)
ax.set_xticks(minor_x_ticks, minor=True)
ax.set_yticks(major_y_ticks)
ax.set_yticks(minor_y_ticks, minor=True)
ax.grid(True, which = 'major', linestyle='-', color='black')
ax.grid(True, which = 'minor', linestyle='--')

plt.xticks(fontsize = tick_fontsize)
plt.yticks(fontsize = tick_fontsize)
plt.xlabel('Freq (GHz)', fontsize=axes_label_fontsize)
plt.ylabel('Impedance', fontsize=axes_label_fontsize)
plt.title("Impedance vs Freq (GHz)", fontsize = title_fontsize)
plt.plot(freq_GHz, Z_Im_antenna, label='Im(Z)')
plt.plot(freq_GHz, Z_Re_antenna, label='Re(Z)')
plt.legend(loc="center left", fontsize=legend_fontsize)
#plt.savefig("Impedance.csv")