#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 21:57:06 2021

@author: akhil
"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import numpy as np
from io import StringIO
import subprocess
import sys
import time
import os
import shutil
import pandas as pd

r=0.9
fs = 24.0
s_size=4.0
s_space=1
lwidth=1.5

preamble='\\usepackage{times}\n\\usepackage{newtxmath}\n\\usepackage{siunitx}\n'
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=preamble)
# Say, "the default sans-serif font is COMIC SANS"
matplotlib.rcParams['font.serif'] = "Times New Roman"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "serif"



C = np.genfromtxt('Dump_Mesogens/Director_11.xyz',delimiter=',',skip_header=2)

fig, ax = plt.subplots(figsize=(7,7),dpi=100,subplot_kw={'projection':'polar'})

# ax.hist(Angles,18,weights=C[j][4:22])
ax.scatter(C[:,18]*np.pi/180.0, C[:,19],s=C[:,19]**2,c=C[:,18]*np.pi/180.0,cmap='hsv',alpha=0.5)  # Plot some data on the axes.
# ax.plot(Dir[:,0], Dir[:,2],'b.',markersize=s_size,markevery=s_space) 
# ax.plot(Dir[:,0], Dir[:,3],'g.',markersize=s_size,markevery=s_space) 

ordertext = 'Order = ' + str(round(np.mean(C[:,15])*100.0)/100.0)
ax.text(66*np.pi/180,23.5,ordertext,color='k',fontsize=fs)

for i in range(0,H.shape[0]):
	ax.text((Angles[i])*np.pi/180.0,18,str(H[i]),color='r',fontsize=fs)


# ax.set_xlim(0.0,5)

# ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.set_rmax(22)
ax.set_rticks([0,5,10,15,20])

ax.set_thetamin(0)
ax.set_thetamax(90)

#Choose X and Y Labels
# ax[0].set_title("a",color= 'k', fontsize= fs,loc='left')
# ax.set_xlabel(r'Time (ns)',color= 'k', fontsize= fs)
# ax.set_ylabel(r'Order parameter', color= 'k', fontsize= fs)

ax.minorticks_on()

for tick in ax.get_xticklabels():
    tick.set_fontsize(r*fs)
for tick in ax.get_yticklabels():
    tick.set_fontsize(r*fs)

ax.tick_params(which='major',direction='in', length=10, width=1.5, colors='k')
ax.tick_params(which='minor',direction='in', length=5, width=1.5, colors='k')
ax.tick_params(which='both',top=True, right=True)

#ax.legend(loc='upper left', bbox_to_anchor=(0.1,0.3), shadow=False, fontsize=r*fs, frameon=False,ncol=1,columnspacing=0.8)
           
ratio=1.0
# ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
#plt.show()
# ax.set_aspect('equal')
plt.savefig(figname,dpi=300);
