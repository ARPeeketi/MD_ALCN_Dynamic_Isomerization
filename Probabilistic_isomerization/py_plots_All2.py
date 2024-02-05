#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 22 Feb 2022

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
lwidth = 2.0
stdalpha = 0.25
wi = 10000
Slevel = 2.0

foldername = ''

print(os.getcwd()[-40:])
filename = foldername+'Data_TTPDLO.dat'
A=np.genfromtxt(filename,skip_header=0,skip_footer=0)

filename = foldername+'Data_Diheds.dat'
Di=np.genfromtxt(filename,skip_header=0,skip_footer=0)
Di[:,1]=Di[:,1]/1000.0
Na = A.shape[0]
# print(Na)
B=np.zeros([Na,1])
for i in range(0,Na):
	B[i,0] = A[i,1]/1000.0 #to convert fs to ps

# print(B[0,0],B[-1,0])
N = Na-1
Ttime = B[-1,0]
wi = 10*10000 #10 ps
Av_density = 1.1206
Av_order = np.mean(A[0:wi,11])

Av_density1 = np.mean(A[-wi:,4])
Av_order1 = np.mean(A[-wi:,11])

Dchange = 100.0*(Av_density-Av_density1)/Av_density
Ochange = 100.0*(Av_order-Av_order1)/Av_order



print('Total Time', ' %8.4f ' % Ttime, 'ps    Total Iters', N)
print('Time        - %12.4f %12.4f' % (np.mean(B[0:wi,0]),np.mean(B[-wi:,0])))
print('Temperature - %12.4f %12.4f' % (np.mean(A[0:wi,2]),np.mean(A[-wi:,2])))
print('Pressure    - %12.4f %12.4f' % (np.mean(A[0:wi,3]),np.mean(A[-wi:,3])))
print('Density     - %12.4f %12.4f' % (np.mean(A[0:wi,4]),np.mean(A[-wi:,4])))
print('%% Density change  - %12.4f' % (Dchange))
print('Order       - %12.4f %12.4f' % (np.mean(A[0:wi,11]),np.mean(A[-wi:,11])))
print('%% Order change    - %12.4f' % (Ochange))
print('Avg Dihed   - %12.4f %12.4f' % (np.mean(Di[0:wi,12]),np.mean(Di[-wi:,12])))
print('No. of cis  - %12.4f %12.4f' % (np.mean(Di[0:wi,14]),np.mean(Di[-wi:,14])))
print('Dimension_X - %12.4f %12.4f' % (np.mean(A[0:wi,8]),np.mean(A[-wi:,8])))
print('Dimension_Y - %12.4f %12.4f' % (np.mean(A[0:wi,9]),np.mean(A[-wi:,9])))
print('Dimension_Z - %12.4f %12.4f' % (np.mean(A[0:wi,10]),np.mean(A[-wi:,10])))


# print('\n')
# for j in range(2,12):
	# print('Dihed angle - %12.4f %12.4f' % (np.max(Di[100000:,j]),np.min(Di[100000:,j]) ))

foldername = os.getcwd()
foldername2=foldername.split('/')
homedire = '/home/Akhil_CRSM/Akhil/LAMMPS/Oct_2022/M10D1A1_10X/Final/'
filename = homedire + 'All_data/' + foldername2[-2] + '_' + foldername2[-1]+ '.dat'
fid = open(filename,'w')
line = 'Total Time %8.4f ps' % Ttime +  ' Total Iters %d\n' % N
fid.write(line)
line = 'Time        - %12.4f %12.4f\n' % (np.mean(B[:wi,0]),np.mean(B[-wi-1:,0]))
fid.write(line)
line = 'Temperature - %12.4f %12.4f\n' % (np.mean(A[:wi,2]),np.mean(A[-wi-1:,2]))
fid.write(line)
line = 'Pressure    - %12.4f %12.4f\n' % (np.mean(A[:wi,3]),np.mean(A[-wi-1:,3]))
fid.write(line)
line = 'Density     - %12.4f %12.4f\n' % (np.mean(A[:wi,4]),np.mean(A[-wi-1:,4]))
fid.write(line)
# line = '%% Density change  - %12.4f\n' % (Dchange)
# fid.write(line)
line = 'Order       - %12.4f %12.4f\n' % (np.mean(A[:wi,11]),np.mean(A[-wi-1:,11]))
fid.write(line)
# line = '%% Order change    - %12.4f\n' % (Ochange)
# fid.write(line)
line = 'Avg Dihed   - %12.4f %12.4f\n' % (np.mean(Di[:wi,12]),np.mean(Di[-wi-1:,12]))
fid.write(line)
line = 'No. of cis  - %12.4f %12.4f\n' % (np.mean(Di[:wi,14]),np.mean(Di[-wi-1:,14]))
fid.write(line)
line = 'Dimension_X - %12.4f %12.4f\n' % (np.mean(A[:wi,8]),np.mean(A[-wi-1:,8]))
fid.write(line)
line = 'Dimension_Y - %12.4f %12.4f\n' % (np.mean(A[:wi,9]),np.mean(A[-wi-1:,9]))
fid.write(line)
line = 'Dimension_Z - %12.4f %12.4f\n' % (np.mean(A[:wi,10]),np.mean(A[-wi-1:,10]))
fid.write(line)
fid.close()

foldername = ''
#------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Fig Temperature
wi = 20000
figname= foldername+'Fig_Temperature.PNG'
fig, ax = plt.subplots(figsize=(7.0,7),dpi=50)
#Choose range
Amin = np.min(A[:,2])
Amax = np.max(A[:,2])

Amin = np.round(Amin - abs(Amax-Amin)/10,0)
Amax = np.round(Amax + abs(Amax-Amin)/10,0)


Ddata = pd.Series(A[:,2])


MV = Ddata.rolling(wi).mean()
STD = Ddata.rolling(wi).std()

MVa = MV.to_numpy()
STDa = STD.to_numpy()

St = MVa + Slevel*STDa
Sb = MVa - Slevel*STDa

# ax.plot(B[:,0], A[:,2],'r.',markersize=s_size,markevery=s_space)  # Plot some data on the axes.
ax.fill_between(B[:,0], St,Sb,facecolor='red',alpha=stdalpha)  # Plot some data on the axes.
ax.plot(B[:,0], MVa,'r-',markersize=s_size,markevery=s_space, linewidth=lwidth)  # Plot some data on the axes.


ax.set_xlim(0.0,Ttime)

# ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
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
figname= foldername+'Fig_Pressure.PNG'
fig, ax = plt.subplots(figsize=(7,7),dpi=50)
#Choose range
Amin = np.min(A[:,3])
Amax = np.max(A[:,3])

Amin = np.round(Amin - abs(Amax-Amin)/10,0)
Amax = np.round(Amax + abs(Amax-Amin)/10,0)

Ddata = pd.Series(A[:,3])

MV = Ddata.rolling(wi).mean()
STD = Ddata.rolling(wi).std()

MVa = MV.to_numpy()
STDa = STD.to_numpy()

St = MVa + Slevel*STDa
Sb = MVa - Slevel*STDa
ax.fill_between(B[:,0], St,Sb,facecolor='blue',alpha=stdalpha)  # Plot some data on the axes.
ax.plot(B[:,0], MVa,'b-',markersize=s_size,markevery=s_space, linewidth=lwidth)  # Plot some data on the axes.


ax.set_xlim(0.0,Ttime)

# ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.set_ylim(Amin,Amax)

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
figname= foldername+'Fig_Density.PNG'
fig, ax = plt.subplots(figsize=(7.0,7),dpi=50)
#Choose range
Amin = np.min(A[:,4])
Amax = np.max(A[:,4])

Amin = np.round(10000*(Amin - abs(Amax-Amin)/10),0)/10000
Amax = np.round(10000*(Amax + abs(Amax-Amin)/10),0)/10000

Ddata = pd.Series(A[:,4])


MV = Ddata.rolling(wi).mean()
STD = Ddata.rolling(wi).std()

MVa = MV.to_numpy()
STDa = STD.to_numpy()

St = MVa + Slevel*STDa
Sb = MVa - Slevel*STDa
ax.fill_between(B[:,0], St,Sb,facecolor='red',alpha=stdalpha)  # Plot some data on the axes.
ax.plot(B[:,0], MVa,'r-',markersize=s_size,markevery=s_space, linewidth=lwidth)  # Plot some data on the axes.


ax.set_xlim(0.0,Ttime)
# ax.set_xticks(np.linspace())
# ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.set_ylim(Amin,Amax)


#Choose X and Y Labels
# ax[0].set_title("a",color= 'k', fontsize= fs,loc='left')
ax.set_xlabel(r'Time (\SI{}{\pico\second})',color= 'k', fontsize= fs)
ax.set_ylabel(r'Density (\SI{}{\gram\per\cubic\centi\meter})', color= 'r', fontsize= fs)

ax.minorticks_on()

for tick in ax.get_xticklabels():
    tick.set_fontsize(r*fs)
for tick in ax.get_yticklabels():
    tick.set_fontsize(r*fs)

ax.tick_params(which='major',axis='x',direction='in', length=10, width=1.5, colors='k')
ax.tick_params(which='minor',axis='x',direction='in', length=5, width=1.5, colors='k')
ax.tick_params(which='both',top=True, right=False)

ax.tick_params(which='major',axis='y',direction='in', length=10, width=1.5, colors='r')
ax.tick_params(which='minor',axis='y',direction='in', length=5, width=1.5, colors='r')

ax.spines['left'].set_color('red')
ax2=ax.twinx()
C=B[wi-1:,0]
MVa=MVa[wi-1:]
MVaP = (Av_density-MVa[:])*100.0/Av_density

Amin = np.min(MVaP)
Amax = np.max(MVaP)
Amin = np.round(Amin - abs(Amax-Amin)/10,0)
Amax = np.round(Amax + abs(Amax-Amin)/10,0)

ax2.plot(C, MVaP,'b-',markersize=s_size,markevery=s_space, linewidth=lwidth) 
ax2.set_ylim(Amin,Amax)
ax2.set_ylabel(r'\% Decrease in density', color= 'b', fontsize= fs)
ax2.tick_params(which='major',axis='y',direction='in', length=10, width=1.5, colors='b')
ax2.tick_params(which='minor',axis='y',direction='in', length=5, width=1.5, colors='b')

ax2.minorticks_on()
for tick in ax2.get_yticklabels():
    tick.set_fontsize(r*fs)
ax2.spines['left'].set_color('red')
ax2.spines['right'].set_color('blue')
ratio=1.0
# ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
#plt.show()
# ax.set_aspect('equal')
plt.savefig(figname,bbox_inches='tight',dpi=300);

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Fig Order
figname= foldername+'Fig_Order.PNG'
fig, ax = plt.subplots(figsize=(7.0,7),dpi=50)
#Choose range
Amin = np.min(A[:,11])*100
Amax = np.max(A[:,11])*100

Amin = np.round((Amin - abs(Amax-Amin)/10),0)/100
Amax = np.round((Amax + abs(Amax-Amin)/10),0)/100

Ddata = pd.Series(A[:,11])


MV = Ddata.rolling(wi).mean()
STD = Ddata.rolling(wi).std()

MVa = MV.to_numpy()
STDa = STD.to_numpy()

St = MVa + Slevel*STDa
Sb = MVa - Slevel*STDa
ax.fill_between(B[:,0], St,Sb,facecolor='black',alpha=stdalpha)  # Plot some data on the axes.
ax.plot(B[:,0], MVa,'k-',markersize=s_size,markevery=s_space, linewidth=lwidth)  # Plot some data on the axes.


ax.set_xlim(0.0,Ttime)
# ax.set_xticks(np.linspace())
# ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.set_ylim(Amin,Amax)


#Choose X and Y Labels
# ax[0].set_title("a",color= 'k', fontsize= fs,loc='left')
ax.set_xlabel(r'Time (\SI{}{\pico\second})',color= 'k', fontsize= fs)
ax.set_ylabel(r'Order', color= 'r', fontsize= fs)

ax.minorticks_on()

for tick in ax.get_xticklabels():
    tick.set_fontsize(r*fs)
for tick in ax.get_yticklabels():
    tick.set_fontsize(r*fs)

ax.tick_params(which='major',axis='x',direction='in', length=10, width=1.5, colors='k')
ax.tick_params(which='minor',axis='x',direction='in', length=5, width=1.5, colors='k')
ax.tick_params(which='both',top=True, right=False)

ax.tick_params(which='major',axis='y',direction='in', length=10, width=1.5, colors='r')
ax.tick_params(which='minor',axis='y',direction='in', length=5, width=1.5, colors='r')

ax.spines['left'].set_color('red')
ax2=ax.twinx()
C=B[wi-1:,0]
MVa=MVa[wi-1:]
MVaP = (Av_order-MVa[:])*100.0/Av_order

Amin = np.min(MVaP)
Amax = np.max(MVaP)
Amin = np.round(Amin - abs(Amax-Amin)/10,0)
Amax = np.round(Amax + abs(Amax-Amin)/10,0)

ax2.plot(C, MVaP,'b-',markersize=s_size,markevery=s_space, linewidth=lwidth) 
ax2.set_ylim(Amin,Amax)
ax2.set_ylabel(r'\% Change in order', color= 'b', fontsize= fs)
ax2.tick_params(which='major',axis='y',direction='in', length=10, width=1.5, colors='b')
ax2.tick_params(which='minor',axis='y',direction='in', length=5, width=1.5, colors='b')

ax2.minorticks_on()
for tick in ax2.get_yticklabels():
    tick.set_fontsize(r*fs)
ax2.spines['left'].set_color('red')
ax2.spines['right'].set_color('blue')
ratio=1.0
# ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
#plt.show()
# ax.set_aspect('equal')
plt.savefig(figname,bbox_inches='tight',dpi=300);


#------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Fig Avg Dihed angle
figname= foldername+'Fig_Avg_Dihed.PNG'
fig, ax = plt.subplots(figsize=(7.0,7),dpi=50)
#Choose range
Amin = np.min(Di[:,12])
Amax = np.max(Di[:,12])

Amin = np.round(Amin - abs(Amax-Amin)/10,0)
Amax = np.round(Amax + abs(Amax-Amin)/10,0)


Ddata = pd.Series(Di[:,12])


MV = Ddata.rolling(wi).mean()
STD = Ddata.rolling(wi).std()

MVa = MV.to_numpy()
STDa = STD.to_numpy()

St = MVa + Slevel*STDa
Sb = MVa - Slevel*STDa

# ax.plot(B[:,0], A[:,2],'r.',markersize=s_size,markevery=s_space)  # Plot some data on the axes.
ax.fill_between(Di[:,1], St,Sb,facecolor='red',alpha=stdalpha)  # Plot some data on the axes.
ax.plot(Di[:,1], MVa,'r-',markersize=s_size,markevery=s_space, linewidth=lwidth)  # Plot some data on the axes.


ax.set_xlim(0.0,Ttime)

# ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.set_ylim(0,180.0)
ax.set_yticks(np.arange(0,180.1,30))
#Choose X and Y Labels
# ax[0].set_title("a",color= 'k', fontsize= fs,loc='left')
ax.set_xlabel(r'Time (\SI{}{\pico\second})',color= 'k', fontsize= fs)
ax.set_ylabel(r'Avg. Dihed Angle (\SI{}{\degree})', color= 'k', fontsize= fs)

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
#----------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
# Fig Temperature
figname=foldername+ 'Fig_Cis_MF.PNG'
fig, ax = plt.subplots(figsize=(7.0,7),dpi=50)
#Choose range

Amin = np.min(Di[:,14])
Amax = np.max(Di[:,14])

Amin = np.round(Amin - abs(Amax-Amin)/10,0)
Amax = np.round(Amax + abs(Amax-Amin)/10,0)


Ddata = pd.Series(Di[:,14])


MV = Ddata.rolling(wi).mean()
STD = Ddata.rolling(wi).std()

MVa = MV.to_numpy()
STDa = STD.to_numpy()

St = MVa + Slevel*STDa
Sb = MVa - Slevel*STDa

# ax.plot(B[:,0], A[:,2],'r.',markersize=s_size,markevery=s_space)  # Plot some data on the axes.
ax.fill_between(Di[:,1], St,Sb,facecolor='red',alpha=stdalpha)  # Plot some data on the axes.
ax.plot(Di[:,1], MVa,'r-',markersize=s_size,markevery=s_space, linewidth=lwidth)  # Plot some data on the axes.


ax.set_xlim(0.0,Ttime)

# ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.set_ylim(-1,11)
ax.set_yticks(np.arange(0,10.1,2))
#Choose X and Y Labels
# ax[0].set_title("a",color= 'k', fontsize= fs,loc='left')
ax.set_xlabel(r'Time (\SI{}{\pico\second})',color= 'k', fontsize= fs)
ax.set_ylabel(r'No. of \emph{cis} isomers', color= 'k', fontsize= fs)

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
# Fig Box Dimensions
figname= foldername+'Fig_Length_XYZ.PNG'
fig, ax = plt.subplots(figsize=(7,7),dpi=50)
#Choose range
Amin = np.min(np.min(A[:,8:11]))
Amax = np.max(np.max(A[:,8:11]))

Amin = np.round(10000*(Amin - abs(Amax-Amin)/10),0)/10000
Amax = 1.01*np.round(10000*(Amax + abs(Amax-Amin)/10),0)/10000


Ddata = pd.Series(A[:,8])

MV = Ddata.rolling(wi).mean()
STD = Ddata.rolling(wi).std()

MVa = MV.to_numpy()
STDa = STD.to_numpy()

St = MVa + Slevel*STDa
Sb = MVa - Slevel*STDa
ax.fill_between(B[:,0], St,Sb,facecolor='red',alpha=stdalpha)  # Plot some data on the axes.
ax.plot(B[:,0], MVa,'r-',markersize=s_size,markevery=s_space, linewidth=lwidth,label='X')  # Plot some data on the axes.

Ddata = pd.Series(A[:,9])


MV = Ddata.rolling(wi).mean()
STD = Ddata.rolling(wi).std()

MVa = MV.to_numpy()
STDa = STD.to_numpy()

St = MVa + Slevel*STDa
Sb = MVa - Slevel*STDa
ax.fill_between(B[:,0], St,Sb,facecolor='blue',alpha=stdalpha)  # Plot some data on the axes.
ax.plot(B[:,0], MVa,'b-',markersize=s_size,markevery=s_space, linewidth=lwidth,label='Y')  # Plot some data on the axes.

Ddata = pd.Series(A[:,10])

MV = Ddata.rolling(wi).mean()
STD = Ddata.rolling(wi).std()

MVa = MV.to_numpy()
STDa = STD.to_numpy()

St = MVa + Slevel*STDa
Sb = MVa - Slevel*STDa
ax.fill_between(B[:,0], St,Sb,facecolor='green',alpha=stdalpha)  # Plot some data on the axes.
ax.plot(B[:,0], MVa,'g-',markersize=s_size,markevery=s_space, linewidth=lwidth,label='Z')  # Plot some data on the axes.


ax.set_xlim(0.0,Ttime)

# ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

ax.set_ylim(Amin,Amax)

#Choose X and Y Labels
# ax[0].set_title("a",color= 'k', fontsize= fs,loc='left')
ax.set_xlabel(r'Time (\SI{}{\pico\second})',color= 'k', fontsize= fs)
ax.set_ylabel(r'Unit cell dimensions ($\AA$)', color= 'k', fontsize= fs)

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
