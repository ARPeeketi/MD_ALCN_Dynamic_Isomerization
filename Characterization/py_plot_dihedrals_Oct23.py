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
import pandas as pd
import shutil
import os

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
lwidth = 1.0
stdalpha = 0.25
wi = 1000
Slevel = 1.0

filename = 'Data_Diheds.dat'
Di=np.genfromtxt(filename,skip_header=0,skip_footer=0)
Di[:,1]=Di[:,1]/1000.0 - 10.0

N = Di.shape[0]
A = np.zeros((N,11))
B = np.zeros((N,11))

for i in range(0,10):
	Data = pd.Series(np.absolute(Di[:,i+2]))
	A[:,i] = Data.rolling(wi).mean()
	B[:,i] = Data.rolling(wi).std()

for i in range(1,11):
	figname= 'Diheds/Fig_Dihed' + str(i) + '.PNG'
	fig, ax = plt.subplots(figsize=(7.0,7.0),dpi=50)


	ax.plot(Di[:,1], A[:,i-1],'r',linewidth=lwidth)  # Plot some data on the axes.

	ax.set_xlim(-10.0,100.0)
	ax.set_xticks(np.array([0,20,40,60,80,100]))
	# ~ ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
	ax.set_ylim(-3,183.0)
	ax.set_yticks(np.arange(0,180.1,30))
	#Choose X and Y Labels
	# ax[0].set_title("a",color= 'k', fontsize= fs,loc='left')
	ax.set_xlabel(r'Time (\SI{}{\pico\second})',color= 'k', fontsize= fs)
	ax.set_ylabel(r'$\phi$ (\SI{}{\degree})', color= 'k', fontsize= fs)

	ax.minorticks_on()

	for tick in ax.get_xticklabels():
		tick.set_fontsize(r*fs)
	for tick in ax.get_yticklabels():
		tick.set_fontsize(r*fs)

	ax.tick_params(which='major',direction='in', length=10, width=1.5, colors='k')
	ax.tick_params(which='minor',direction='in', length=5, width=1.5, colors='k')
	ax.tick_params(which='both',top=True, right=True)

	ratio=1.0
	# ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
	#plt.show()
	# ax.set_aspect('equal')
	plt.savefig(figname,bbox_inches='tight',dpi=300);	
	
for i in range(1,11):
	figname= 'Diheds/Zoom_Fig_Dihed' + str(i) + '.PNG'
	fig, ax = plt.subplots(figsize=(7.0,7.0),dpi=50)


	ax.plot(Di[:,1], A[:,i-1],'r',linewidth=lwidth)  # Plot some data on the axes.

	ax.set_xlim(50.0,55.0)
	ax.set_xticks(np.array([50,51,52,53,54,55]))
	# ~ ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
	ax.set_ylim(-3,183.0)
	ax.set_yticks(np.arange(0,180.1,30))
	#Choose X and Y Labels
	# ax[0].set_title("a",color= 'k', fontsize= fs,loc='left')
	ax.set_xlabel(r'Time (\SI{}{\pico\second})',color= 'k', fontsize= fs)
	ax.set_ylabel(r'$\phi$ (\SI{}{\degree})', color= 'k', fontsize= fs)

	ax.minorticks_on()

	for tick in ax.get_xticklabels():
		tick.set_fontsize(r*fs)
	for tick in ax.get_yticklabels():
		tick.set_fontsize(r*fs)

	ax.tick_params(which='major',direction='in', length=10, width=1.5, colors='k')
	ax.tick_params(which='minor',direction='in', length=5, width=1.5, colors='k')
	ax.tick_params(which='both',top=True, right=True)

	ratio=1.0
	# ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
	#plt.show()
	# ax.set_aspect('equal')
	plt.savefig(figname,bbox_inches='tight',dpi=300);	
