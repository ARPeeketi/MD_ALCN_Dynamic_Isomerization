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
import imageio

r=0.9
fs = 24.0
s_size=10.0
s_space=1
lwidth = 1.0

preamble='\\usepackage{times}\n\\usepackage{newtxmath}\n\\usepackage{siunitx}\n'
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=preamble)
# Say, "the default sans-serif font is COMIC SANS"
matplotlib.rcParams['font.serif'] = "Times New Roman"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "serif"

C=np.loadtxt('Diheds.dat',delimiter=' ')

# Dihedral angle figure count
figname= 'Dihedrals_t_c.PNG'

fig, ax = plt.subplots(figsize=(7.5,7),dpi=50)

#Choose range
Amin = 0.4
Amax = 0.8

C = C[C[:,0]>0,:]
# ax.plot(C[:,0], C[:,23],'r-',markersize=s_size,markevery=s_space,linewidth=lwidth,label='\emph{trans}') 
ax.plot(C[:,0]-200, C[:,22]/10.0,'b-',markersize=s_size,markevery=s_space,linewidth=lwidth,label = '\emph{cis}')  

ax.set_xlim(0.0,600.0)
ax.set_xticks(np.arange(0.0,600.1,100.0))
# ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.set_ylim(-0.10,1.10)
ax.set_yticks(np.arange(0.0,1.01,0.2))

#Choose X and Y Labels
# ax[0].set_title("a",color= 'k', fontsize= fs,loc='left')
ax.set_xlabel(r'Time (\SI{}{\pico\second})',color= 'k', fontsize= fs)
ax.set_ylabel(r'$n_{c}$', color= 'k', fontsize= fs)

ax.minorticks_on()

for tick in ax.get_xticklabels():
    tick.set_fontsize(r*fs)
for tick in ax.get_yticklabels():
    tick.set_fontsize(r*fs)
    tick.set_fontsize(r*fs)

ax.tick_params(which='major',direction='in', length=10, width=1.5, colors='k')
ax.tick_params(which='minor',direction='in', length=5, width=1.5, colors='k')
ax.tick_params(which='both',top=True, right=True)

# ax.legend(loc='upper left', bbox_to_anchor=(0.1,1.01), shadow=False, fontsize=r*fs, frameon=False,ncol=2,columnspacing=0.8)
           
ratio=1.0
ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
# plt.show()
# ax.set_aspect('equal')
plt.savefig(figname,bbox_inches='tight',dpi=300);
