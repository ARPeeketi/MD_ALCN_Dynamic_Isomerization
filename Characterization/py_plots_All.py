#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 21:57:06 2021

@author: akhil
"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import numpy as np
from io import StringIO
import imageio
import pandas as pd
import shutil
import os

from scipy.interpolate import interp1d, make_interp_spline, BSpline
from scipy.fft import fft, fftfreq
from scipy.signal import savgol_filter

preamble='\\usepackage{times}\n\\usepackage{newtxmath}\n\\usepackage{siunitx}\n'
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=preamble)
# Say, "the default sans-serif font is COMIC SANS"
matplotlib.rcParams['font.serif'] = "Times New Roman"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "serif"

r=0.9
fs = 24.0
s_size=4.0
s_space=1
lwidth = 2.0
stdalpha = 0.25
wi = 10000
Slevel = 2.0




filename = 'Data_C200ps.dat'
A=np.genfromtxt(filename,delimiter=',')


Na = A.shape[0]
B = np.zeros([Na,1])
print(Na)
for i in range(0,Na):
	B[i,0] = i*1.0*0.10*1.0e-3

print(B[0,0],B[-1,0])

# print(B[0,0],B[-1,0])
N = Na-1
Ttime = B[-1,0]

wi = 10000


#----------------------------------------------------------------------------------------------------------------
# Fig Density
figname= 'Fig_Density2.PNG'
fig, ax = plt.subplots(figsize=(7.0,7),dpi=50)
#Choose range
Amin = np.min(A[:,4])
Amax = np.max(A[:,4])

Amin = np.round(10000*(Amin - abs(Amax-Amin)/10),0)/10000
Amax = np.round(10000*(Amax + abs(Amax-Amin)/10),0)/10000
po=3
wl = 20001
y_smooth = savgol_filter(A[:,4], window_length=wl, polyorder=po, mode="interp")


A1 = np.hstack(([-0.00001],B[:,0]))
A2 = np.hstack(([A[0,4]],y_smooth))
# ax.plot([100,100],[0.8,1.2],'k--',linewidth=1.0)
ax.plot(A1, A2,'r-',markersize=s_size,markevery=s_space, linewidth=lwidth)  # Plot some data on the axes.

# Ddata = pd.Series(A[:,4])


# MV = Ddata.rolling(wi).mean()
# STD = Ddata.rolling(wi).std()

# MVa = MV.to_numpy()
# STDa = STD.to_numpy()

# St = MVa + Slevel*STDa
# Sb = MVa - Slevel*STDa
# ax.fill_between(B[:,0]-200, St,Sb,facecolor='red',alpha=stdalpha)  # Plot some data on the axes.
# ax.plot(B[:,0]-200, MVa,'r-',markersize=s_size,markevery=s_space, linewidth=lwidth)  # Plot some data on the axes.
ax2xlim=0
ax.add_patch(Rectangle((ax2xlim, 0.8), 200, 0.4, linewidth=0.01, edgecolor='k', facecolor=(0,0,1,0.1)))
ax.add_patch(Rectangle((ax2xlim+200, 0.8), 200, 0.4, linewidth=0.01, edgecolor='k', facecolor=(1,0,0,0.1)))
ax.add_patch(Rectangle((ax2xlim+400, 0.8),200, 0.4, linewidth=0.01, edgecolor='k', facecolor=(0,0,1,0.1)))
# ax.plot()
ax.text(40,1.0,r'$trans$',fontsize=fs,color='k')
ax.text(250,1.0,r'$cis$',fontsize=fs,color='k')
ax.text(440,1.0,r'$trans$',fontsize=fs,color='k')
ax.set_xlim(-1.0,600)
ax.yaxis.set_ticks(np.arange(0., 600.01, 200.0))
# ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))

ax.set_ylim(0.89,1.16)

ax.yaxis.set_ticks(np.arange(0.9, 1.151, 0.05))


#Choose X and Y Labels
# ax[0].set_title("a",color= 'k', fontsize= fs,loc='left')
ax.set_xlabel(r'Time (\SI{}{\pico\second})',color= 'k', fontsize= fs)
ax.set_ylabel(r'Density (\SI{}{\gram\per\cubic\centi\meter})', color= 'k', fontsize= fs)

ax.minorticks_on()

for tick in ax.get_xticklabels():
    tick.set_fontsize(r*fs)
for tick in ax.get_yticklabels():
    tick.set_fontsize(r*fs)

ax.tick_params(which='major',axis='x',direction='in', length=10, width=1.5, colors='k')
ax.tick_params(which='minor',axis='x',direction='in', length=5, width=1.5, colors='k')
ax.tick_params(which='both',top=True, right=True)

ax.tick_params(which='major',axis='y',direction='in', length=10, width=1.5, colors='k')
ax.tick_params(which='minor',axis='y',direction='in', length=5, width=1.5, colors='k')

ratio=1.0
# ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
#plt.show()
# ax.set_aspect('equal')
plt.savefig(figname,bbox_inches='tight',dpi=300);


