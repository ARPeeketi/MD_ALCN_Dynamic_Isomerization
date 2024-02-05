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

NNN = 5
r=0.9
fs = 24.0
s_size=4.0
s_space=1

preamble='\\usepackage{times}\n\\usepackage{newtxmath}\n\\usepackage{siunitx}\n'
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=preamble)
# Say, "the default sans-serif font is COMIC SANS"
matplotlib.rcParams['font.serif'] = "Times New Roman"
# Then, "ALWAYS use sans-serif fonts"
matplotlib.rcParams['font.family'] = "serif"

i=1
filename = 'Cycles/Cycle_' + str(i) +'/log/l_npt_ex.log'
A=np.genfromtxt(filename,skip_header=26,skip_footer=38)

for i in range(2,NNN):
	print(i)
	filename = 'Cycles/Cycle_' + str(i) +'/log/l_npt_ex.log'
	A1=np.genfromtxt(filename,skip_header=27,skip_footer=38)
	A = np.vstack((A,A1))
# print(A[:,0:3])

file1name = 'Data_C200ps.dat'
f1=open(file1name,'w')
np.savetxt(f1,A,fmt='%15.6f',delimiter=',')
f1.close()

Na = A.shape[0]

B = np.zeros([Na,1])

for i in range(0,Na):
	B[i,0] = i*1.0*0.10*1.0e-3

# print(B)
Timestep = 0.10


#A[:,1] = A[:,1]/1000000 #convert to ns
A5=A
A6=A
N = int(A[-1][0])
Ttime = A[-1,1]

print('Total Time', Ttime, 'ns    Total Iters', N)
print('       Temperature    Pressure      Density')
print('NPT ex %12.4f %12.4f %12.4f' % (np.mean(A5[0:10,2]),np.mean(A5[0:10,3]),np.mean(A5[0:10,4])))
print('NVE ex %12.4f %12.4f %12.4f' % (np.mean(A6[-10:-1,2]),np.mean(A6[-10:-1,3]),np.mean(A6[-10:-1,4])))

#----------------------------------------------------------------------------------------------------------------
# Fig Temperature
figname= 'Temperature.PNG'
fig, ax = plt.subplots(figsize=(6,6),dpi=50)
#Choose range
Amin = np.min(A[:,2])
Amax = np.max(A[:,2])

Amin = np.round(Amin - abs(Amax-Amin)/10,0)
Amax = np.round(Amax + abs(Amax-Amin)/10,0)

ax.plot(B[:,0], A[:,2],'r.',markersize=s_size,markevery=s_space)  # Plot some data on the axes.




#ax.set_xlim(0.0,Ttime)

#ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.set_ylim(Amin,Amax)

#Choose X and Y Labels
# ax[0].set_title("a",color= 'k', fontsize= fs,loc='left')
ax.set_xlabel(r'Time (\SI{}{\pico\second})',color= 'k', fontsize= fs)
ax.set_ylabel(r'Temperature (\SI{}{\kelvin})', color= 'k', fontsize= fs)

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

#----------------------------------------------------------------------------------------------------------------
# Fig Pressure
figname= 'Pressure.PNG'
fig, ax = plt.subplots(figsize=(6,6),dpi=50)
#Choose range
Amin = np.min(A[:,3])
Amax = np.max(A[:,3])

Amin = np.round(Amin - abs(Amax-Amin)/10,0)
Amax = np.round(Amax + abs(Amax-Amin)/10,0)

ax.plot(B[:,0], A[:,3],'r.',markersize=s_size,markevery=s_space)  # Plot some data on the axes.


#ax.set_xlim(0.0,Ttime)

#ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#ax.set_ylim(Amin,Amax)

#Choose X and Y Labels
# ax[0].set_title("a",color= 'k', fontsize= fs,loc='left')
ax.set_xlabel(r'Time (\SI{}{\pico\second})',color= 'k', fontsize= fs)
ax.set_ylabel(r'Pressure (atm)', color= 'k', fontsize= fs)

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

#----------------------------------------------------------------------------------------------------------------
# Fig Density
figname= 'Density.PNG'
fig, ax = plt.subplots(figsize=(6,6),dpi=50)
#Choose range
Amin = np.min(A[:,4])
Amax = np.max(A[:,4])

Amin = np.round(10000*(Amin - abs(Amax-Amin)/10),0)/10000
Amax = np.round(10000*(Amax + abs(Amax-Amin)/10),0)/10000

ax.plot(B[:,0], A[:,4],'r.',markersize=s_size,markevery=s_space)  # Plot some data on the axes.




#ax.set_xlim(0.0,Ttime)

#ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.set_ylim(0.9,1.15)
ax.set_yticks([0.9,0.95,1.00,1.05,1.10,1.15])
#Choose X and Y Labels
# ax[0].set_title("a",color= 'k', fontsize= fs,loc='left')
ax.set_xlabel(r'Time (\SI{}{\pico\second})',color= 'k', fontsize= fs)
ax.set_ylabel(r'Density (\SI{}{\gram\per\cubic\centi\meter})', color= 'k', fontsize= fs)

ax.minorticks_on()

for tick in ax.get_xticklabels():
    tick.set_fontsize(r*fs)
for tick in ax.get_yticklabels():
    tick.set_fontsize(r*fs)

ax.tick_params(which='major',direction='in', length=10, width=1.5, colors='k')
ax.tick_params(which='minor',direction='in', length=5, width=1.5, colors='k')
ax.tick_params(which='both',top=True, right=True)

# ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
#plt.show()
# ax.set_aspect('equal')
plt.savefig(figname,bbox_inches='tight',dpi=300);

#----------------------------------------------------------------------------------------------------------------
# Fig Box Dimensions
figname= 'Length_XYZ.PNG'
fig, ax = plt.subplots(figsize=(6,6),dpi=50)
#Choose range
Amin = np.min(np.min(A[:,8:11]))
Amax = np.max(np.max(A[:,8:11]))

Amin = np.round(10000*(Amin - abs(Amax-Amin)/10),0)/10000
Amax = 1.01*np.round(10000*(Amax + abs(Amax-Amin)/10),0)/10000

ax.plot(B[:,0], A[:,8],'ro',markersize=s_size,markevery=s_space*2,label='X')  # Plot some data on the axes.
ax.plot(B[:,0], A[:,9],'bs',markersize=s_size,markevery=s_space*10,label='Y')  # Plot some data on the axes.
ax.plot(B[:,0], A[:,10],'gd',markersize=s_size*1.2,markevery=s_space*50,label='Z')  # Plot some data on the axes.
#ax.plot([[-1,1],A5[-1,1]],[Amin,Amax],'k-',linewidth=1.0)
#ax.plot([A6[-1,1],A6[-1,1]],[Amin,Amax],'k-',linewidth=1.0)

#ax.set_xlim(0.0,Ttime)

#ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.set_ylim(Amin,Amax)

#Choose X and Y Labels
# ax[0].set_title("a",color= 'k', fontsize= fs,loc='left')
ax.set_xlabel(r'Time (\SI{}{\pico\second})',color= 'k', fontsize= fs)
ax.set_ylabel(r'Length ($\AA$)', color= 'k', fontsize= fs)

ax.minorticks_on()

for tick in ax.get_xticklabels():
    tick.set_fontsize(r*fs)
for tick in ax.get_yticklabels():
    tick.set_fontsize(r*fs)

ax.tick_params(which='major',direction='in', length=10, width=1.5, colors='k')
ax.tick_params(which='minor',direction='in', length=5, width=1.5, colors='k')
ax.tick_params(which='both',top=True, right=True)

ax.legend(loc='upper left', bbox_to_anchor=(0.04,0.98), shadow=False, fontsize=r*fs, frameon=False,ncol=3,columnspacing=0.8)
           
ratio=1.0
# ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
#plt.show()
# ax.set_aspect('equal')
plt.savefig(figname,bbox_inches='tight',dpi=300);
