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

Tcycle = 1.0 #fs Cyclic dynamics time step
Tstep = 1.0 #fs lammps time step

Totaltime = 100.0 #ps #Cyclic dynamic isomerization total time
Ini_time = 10.0 #ps #Initial equilibrium dynamics time

NDiheddump = 10 #lammps exporting diheds every __ iters
Ndump = 100 #lammps exporting diheds every __ iters

Natoms = 7278

N=int((Totaltime+Ini_time)*1000.0/Tstep)

stdalpha = 0.25
wi = 20
Slevel = 2.0

octype = 11.0
o2type = 10.0
n2type = 8.0



print(os.getcwd()[-40:])
ofolder= os.getcwd()
ofolder = ofolder + '/'

foldername = ofolder + 'Dump_Mesogens'
if os.path.isdir(foldername):
	shutil.rmtree(foldername,ignore_errors=True)
os.mkdir(foldername)

# Fig Director
print('Extracting mesogen data from dump files')

dumpname = 'Dump_All_cycles.dump'

dumpfile = open(dumpname,'r')

count = 0
Nlines = 10

while (count < Nlines):
	count += 1
	line = dumpfile.readline()
	if (count==1):
		line1 = line
	if (count==5):
		line5 = line
	if (count==9):
		line9 = line	

#print(line1,line9)
dumpfile.close()


dumpfile = open(dumpname,'r')
count = 0
count_iter = 0

Niters =  int(N/Ndump)+1
Nlines = Niters*(8)
A = np.zeros((Niters,Natoms,6))
LBox = np.zeros((Niters,9))

Dir = np.zeros((Niters,20))
while (count < Nlines):
	count += 1
	line = dumpfile.readline()
	if (line == line1):
		count_iter +=1
		line_tstep = dumpfile.readline()
		Tstepc = int(line_tstep)
		Dir[count_iter-1,0] = Tstepc
		#print(line,count,count_iter,Tstep)	
	elif (line==line5):
		line = dumpfile.readline()
		xlo=float(line.split(' ')[0])
		xhi=float(line.split(' ')[1])
		
		line = dumpfile.readline()
		ylo=float(line.split(' ')[0])
		yhi=float(line.split(' ')[1])		
		
		line = dumpfile.readline()
		zlo=float(line.split(' ')[0])
		zhi=float(line.split(' ')[1])		
		
		LBox[count_iter-1,0] = xhi - xlo
		LBox[count_iter-1,1] = yhi - ylo
		LBox[count_iter-1,2] = zhi - zlo
		LBox[count_iter-1,3] = xlo
		LBox[count_iter-1,4] = ylo
		LBox[count_iter-1,5] = zlo
		LBox[count_iter-1,6] = xhi
		LBox[count_iter-1,7] = yhi
		LBox[count_iter-1,8] = zhi
	elif (line==line9):
		#print(line,count,count_iter,Tstep)
		for i in range(0,Natoms):
			line = dumpfile.readline()
			line = line[:-1]
			#print(line.split(' '))
			#print(line.split(' '))
			A[count_iter-1,i,0]=int(line.split(' ')[0])
			A[count_iter-1,i,1]=int(line.split(' ')[1])
			A[count_iter-1,i,2]=float(line.split(' ')[2])
			A[count_iter-1,i,3]=float(line.split(' ')[3])
			A[count_iter-1,i,4]=float(line.split(' ')[4])
			A[count_iter-1,i,5]=float(line.split(' ')[5])
			#print(A[i,:])
			#print(A[i,:])
		# print(A[count_iter-1,0,:])
		#print(max(A[count_iter-1,:,4]))
	#print(line)
	#print()
dumpfile.close()	

print('Files read')
print(N,Niters)
# print(A[0,0,:])
for i in range(0,Niters):	
	#print('Total Atoms Min and Max')
	#print(min(A[i,:,4]))
	#print(max(A[i,:,4]))
	
	Aa = A[i,:,:]
	ci = 0
	Aa = Aa[Aa[:,ci].argsort()]	
	B = Aa[Aa[:,1]==octype,:]
	# Aa = Aa[Aa[:,ci].argsort()]
	#B[:,2:5] = B[:,2:5]
	
	# print(A1)
	B = B[B[:,ci].argsort()]
	
	#print(B[:,[0,1,5]])
	Ndir = int(B.shape[0]/2)
	print(i,Ndir)
	#print(Ndir)
	C = np.zeros((Ndir,22))
	E = np.zeros([Ndir,15])
	E1=np.zeros([1,15])
	# print(Ndir)
	# print(LBox[i,:])
	for j in range(0,Ndir):
		# print(j+1)
		C[j,0] = B[2*j,0]
		C[j,1] = B[2*j+1,0]
		# print(Aa[int(C[j,0])-1:int(C[j,1]),:])
		C[j,2:5] = B[2*j+1,2:5]
		C[j,5:8] = B[2*j,2:5]
		# print(C[j,2:8])
		# print(C[j,5:8])
		# print('\n')
		C[j,8:11] = (C[j,2:5]-C[j,5:8])		
		for k in range(8,11):
			CMDa1 = abs(0.50 - C[j,k-6])
			CMDa2 = abs(0.50 - C[j,k-3])
			if (C[j,k] > 0.5):
				if (CMDa1>CMDa2):
					C[j,k-6] = C[j,k-6] - 1.0 #1st atom
				else:
					C[j,k-3] = 1.0 + C[j,k-3] #2nd atom
			elif (C[j,k] < -0.5):
				if (CMDa1>CMDa2):
					C[j,k-6] = C[j,k-6] + 1.0 #1st atom
				else:
					C[j,k-3] = C[j,k-3] - 1.0 #2nd atom			
		C[j,8:11] = (C[j,2:5]-C[j,5:8])
		for k in range(8,11):
			C[j,k] = LBox[i,k-8]*C[j,k] 
		for k in range(2,5):
			C[j,k] = LBox[i,k-2]*C[j,k] + LBox[i,k+1]
		for k in range(5,8):
			C[j,k] = LBox[i,k-5]*C[j,k] + LBox[i,k-2]
		# print(C[j,2:5])
		# print(C[j,5:8])
				
		C[j,11] = np.sqrt(np.sum(np.multiply(C[j,8:11],C[j,8:11])))
		
		# print('dir length',C[j,11])	
		C[j,12] = C[j,8]/C[j,11]
		C[j,13] = C[j,9]/C[j,11]
		C[j,14] = C[j,10]/C[j,11]
		
		# print(C[j,12:15],'\n')
		# if (C[j,12]<0.0):
			# C[j,12:15] = -C[j,12:15]
			
		C[j,15] = (3.0*C[j,12]**2-1.0)/2.0
		C[j,16] = (3.0*C[j,13]**2-1.0)/2.0
		C[j,17] = (3.0*C[j,14]**2-1.0)/2.0
		
		C[j,18] = np.arccos(C[j,12])*180.0/np.pi
		C[j,19] = C[j,11]
		
		C[j,20] = np.argmax(np.absolute([C[j,12:15]]))
	
		
		if (C[j,0]<C[j,1]):
			Adum = Aa[int(C[j,0])-1:int(C[j,1]),:]
		else:
			Adum = Aa[int(C[j,1])-1:int(C[j,0]),:]
		
		count10 = np.shape(Adum[Adum[:,1]==o2type,:])[0]
		count8 = np.shape(Adum[Adum[:,1]==n2type,:])[0]
		# print(Adum[:,0:2],count10,count8)
		if (count10==1.0):
			C[j,21] = 1.0
		elif (count10==2.0):
			C[j,21] = 2.0
		if (count8==2.0):
			# print(j)
			C[j,21]=3.0

	D = C[C[:,21]<=2,:]
	Ndir2=D.shape[0]
	
	a11=0.0
	a22=0.0
	a33=0.0
	a12=0.0
	a23=0.0
	a13=0.0
	
	for j in range(0,Ndir):
		for k in range(0,Ndir):	
			a11=a11 + (3.0*C[j,12]*C[k,12] - 1.0)/2.0
			a22=a22 + (3.0*C[j,13]*C[k,13] - 1.0)/2.0
			a33=a33 + (3.0*C[j,14]*C[k,14] - 1.0)/2.0
			a12=a12 + 3.0*C[j,12]*C[k,13]/2.0
			a23=a23 + 3.0*C[j,13]*C[k,14]/2.0
			a13=a13 + 3.0*C[j,12]*C[k,14]/2.0
	Dir[i,1] = np.mean(C[:,15])
	Dir[i,2] = np.mean(C[:,16])
	Dir[i,3] = np.mean(C[:,17])
	Dir[i,4] = a11/(float(Ndir)**2)
	Dir[i,5] = a22/(float(Ndir)**2)
	Dir[i,6] = a33/(float(Ndir)**2)
	Dir[i,7] = a12/(float(Ndir)**2)
	Dir[i,8] = a23/(float(Ndir)**2)
	Dir[i,9] = a13/(float(Ndir)**2)
	
	
	a11=0.0
	a22=0.0
	a33=0.0
	a12=0.0
	a23=0.0
	a13=0.0
	
	for j in range(0,Ndir2):
		for k in range(0,Ndir2):	
			a11=a11 + (3.0*D[j,12]*D[k,12] - 1.0)/2.0
			a22=a22 + (3.0*D[j,13]*D[k,13] - 1.0)/2.0
			a33=a33 + (3.0*D[j,14]*D[k,14] - 1.0)/2.0
			a12=a12 + 3.0*D[j,12]*D[k,13]/2.0
			a23=a23 + 3.0*D[j,13]*D[k,14]/2.0
			a13=a13 + 3.0*D[j,12]*D[k,14]/2.0
	Dir[i,11] = np.mean(D[:,15])
	Dir[i,12] = np.mean(D[:,16])
	Dir[i,13] = np.mean(D[:,17])
	Dir[i,14] = a11/(float(Ndir)**2)
	Dir[i,15] = a22/(float(Ndir)**2)
	Dir[i,16] = a33/(float(Ndir)**2)
	Dir[i,17] = a12/(float(Ndir)**2)
	Dir[i,18] = a23/(float(Ndir)**2)
	Dir[i,19] = a13/(float(Ndir)**2)
	
	print(i,Dir[i,0],Dir[i,1],Dir[i,11])
	
	


Dir[:,0] = Dir[:,0]*Tstep/1000.0

np.savetxt('Dump_order.dat', Dir, fmt='%15.6f')


foldername = os.getcwd()
foldername2=foldername.split('/')
homedire = '/home/Akhil_CRSM/Akhil/LAMMPS/Oct_2022/M10D1A1_10X/Final/'
filename = homedire + 'Order_data/' + foldername2[-2] + '_' + foldername2[-1]+ '.dat'

f1=open(filename,'w')
np.savetxt(f1,Dir,fmt='%15.6f')
f1.close()

figname= 'Fig_Order3.PNG'

fig, ax = plt.subplots(figsize=(7.5,7),dpi=50)

#Choose range
Amin = 0.4
Amax = 0.8
Ddata = pd.Series(Dir[:,1])

MV = Ddata.rolling(wi).mean()
STD = Ddata.rolling(wi).std()

MVa = MV.to_numpy()
STDa = STD.to_numpy()

St = MVa + Slevel*STDa
Sb = MVa - Slevel*STDa
ax.fill_between(Dir[:,0], St,Sb,facecolor='black',alpha=stdalpha)  # Plot some data on the axes.
ax.plot(Dir[:,0], MVa,'k-',markersize=s_size,markevery=s_space, linewidth=lwidth,label='Y')  # Plot some data on the axes.


Ddata = pd.Series(Dir[:,11])

MV = Ddata.rolling(wi).mean()
STD = Ddata.rolling(wi).std()

MVa = MV.to_numpy()
STDa = STD.to_numpy()

St = MVa + Slevel*STDa
Sb = MVa - Slevel*STDa
ax.fill_between(Dir[:,0], St,Sb,facecolor='red',alpha=stdalpha)  # Plot some data on the axes.
ax.plot(Dir[:,0], MVa,'r-',markersize=s_size,markevery=s_space, linewidth=lwidth,label='Y')  # Plot some data on the axes.

# ax.plot(Dir[:,0], Dir[:,1],'r.',markersize=s_size,markevery=s_space)  # Plot some data on the axes.

ax.set_ylim(0.1,0.6)
ax.set_yticks(np.arange(0.1,0.601,0.1))

#Choose X and Y Labels
# ax[0].set_title("a",color= 'k', fontsize= fs,loc='left')
ax.set_xlabel(r'Time (\SI{}{\pico\second})',color= 'k', fontsize= fs)
ax.set_ylabel(r'Order ($Q$)', color= 'k', fontsize= fs)

ax.minorticks_on()

for tick in ax.get_xticklabels():
    tick.set_fontsize(r*fs)
for tick in ax.get_yticklabels():
    tick.set_fontsize(r*fs)
    tick.set_fontsize(r*fs)

ax.tick_params(which='major',direction='in', length=10, width=1.5, colors='k')
ax.tick_params(which='minor',direction='in', length=5, width=1.5, colors='k')
ax.tick_params(which='both',top=True, right=True)

#ax.legend(loc='upper left', bbox_to_anchor=(0.1,0.3), shadow=False, fontsize=r*fs, frameon=False,ncol=1,columnspacing=0.8)
           
ratio=1.0
# ax.set_aspect(1.0/ax.get_data_ratio()*ratio)
#plt.show()
# ax.set_aspect('equal')
plt.savefig(figname,bbox_inches='tight',dpi=300);
