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
s_size=2.0
s_space=1

preamble='\\usepackage{times}\n\\usepackage{newtxmath}\n\\usepackage{siunitx}\n'
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=preamble)
# Say, "the default sans-serif font is COMIC SANS"
matplotlib.rcParams['font.serif'] = "Times New Roman"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "serif"

Ncycles = 5

j = 0
Angles_1=np.linspace(0,180,19)
Angles=np.linspace(5,180,18)

print('Extracting the dihed hist data')

N = 2000000
NDihed = 1000

Tstep = 0.1
Dstep = NDihed*Tstep/1000.0

NNN = int(Ncycles*N/NDihed)
C = np.zeros([NNN,24])

j=1
for Cycle_num in range(1,5):
	for i in range(NDihed,N+1,NDihed):
		fname =  'Cycles/Cycle_'+ str(Cycle_num) + '/dump/dihed_' + str(i) + '.dat'
		C1=np.genfromtxt(fname,skip_header=9,skip_footer=0)
		C2a=C1[C1[:,4]==81.0,:]
		C2b=C1[C1[:,4]==108.0,:]
		
		if ((len(C2a)>0)and(len(C2b)>0)):
			C2 = np.vstack((C2a,C2b))
		elif ((len(C2a)>0)and(len(C2b)==0)):
			C2 = C2a
		elif ((len(C2a)==0)and(len(C2b)>0)):
			C2 = C2b
		
		ta = j*Dstep
		#C2 = np.vstack([C2a,C2b])
		C2[:,5] = abs(C2[:,5])
		C[j][0] = ta
		C[j][1] = np.mean(C2[:,5])
		C[j][2] = np.min(C2[:,5])
		C[j][3] = np.max(C2[:,5])
		C[j][4:22] = np.histogram([C2[:,5]],bins=Angles_1)[0]
		C[j,22] = np.sum(C[j,4:12])
		C[j,23] = np.sum(C[j,12:22])
		print(ta)
		# fname2 = 'dihed_figures/dihed_' + str(Cycle_num) + '_' + str(int(i/NDihed)) + '.PNG'
		
		# fig, ax = plt.subplots(figsize=(7.0,5.0),dpi=50)
		# ax.hist(Angles,18,weights=C[j][4:22],color='r')
		# ax.set_xlabel(r'Dihedral Angle',color= 'k', fontsize=fs)
		# ax.set_ylabel(r'Frequency', color= 'k', fontsize=fs)
		# ax.set_xlim([0,180])
		# ax.set_ylim([0,10])
		# ax.set_yticks(np.arange(0,10.1,2))
		# ax.minorticks_on()

		# for tick in ax.get_xticklabels():
			# tick.set_fontsize(r*fs)
		# for tick in ax.get_yticklabels():
			# tick.set_fontsize(r*fs)
			# tick.set_fontsize(r*fs)

		# ax.tick_params(which='major',direction='in', length=10, width=1.5, colors='k')
		# ax.tick_params(which='minor',direction='in', length=5, width=1.5, colors='k')
		# ax.tick_params(which='both',top=True, right=True)
		
		# timehere = '%6.2f' % ta
		# titlehere = 'Time ' + timehere + ' ps'
		# plt.title(titlehere, color= 'k', fontsize=fs)
		# plt.savefig(fname2,bbox_inches='tight',dpi=300);
		# plt.close()
		
		#print(C[j][4:22])
		j = j+1

C = C[~np.all(C == 0, axis=1)]
np.savetxt('Diheds.dat', C, fmt=['%0.3f','%.3f','%.3f','%.3f','%d','%d','%d','%d','%d','%d','%d','%d','%d','%d','%d','%d','%d','%d','%d','%d','%d','%d','%d','%d'])
