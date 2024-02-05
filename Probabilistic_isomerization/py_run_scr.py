#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 21:57:06 2021

@author: akhil
"""

import numpy as np
from io import StringIO
import subprocess
import sys
import time
import os
import shutil
import time
from datetime import datetime

starttime = datetime.now()
ofolder= os.getcwd()
ofolder = ofolder + '/'
mdname = 'M10D1A1_C25_eq2'
lammpscmd = '/usr/lib64/openmpi/bin/mpirun -np 8 --hostfile my_hosts /home/Akhil_CRSM/Installations/lammps/src/lmp_mpi -in lammps2.in >> log/terminal.log'
#subprocess.call(lammpscmd, shell=True)

#Dihedrals 81 is thermal, 108 - is tc and 109 is ct
probtc = 0.003
probct = 0.003

Totaltime = 100.0 #ps #Cyclic dynamic isomerization total time
Ini_time = 10.0 #ps #Initial equilibrium dynamics time

Tcycle = 1.0 #fs Cyclic dynamics time step
Tstep = 0.01 #fs lammps time step
N = int(Tcycle/Tstep)
Ncycledelete = 100
Nmesogendump = 100 #Every n number of cycles

NDiheddump = 10 #lammps exporting diheds every __ iters
Ndump = 10 #lammps exporting diheds every __ iters

Ndihedline1 = 35987
Ndiheds = 17334
Ndihedline2 = Ndihedline1 + Ndiheds + 2
TNlines = 59845 +1
Natoms = 7278
Nisomers = 10

octype = 11.0
o2type = 10.0
n2type = 8.0

Actual_iter = 0

dihedic = int(Tcycle/(Tstep*Ndump))



Ncycles1 = int(Ini_time*1000.0/Tcycle) + 1
Dy_Ncycles = int(Totaltime*1000.0/Tcycle)
Ncycles = Dy_Ncycles + Ncycles1

Nlineseach = Natoms + 9
TNlinesdump = Nlineseach*(dihedic+1)

file1name = 'Data_TTPDLO.dat'
file2name = 'Data_Diheds.dat'
dumpfilea = 'Dump_All_cycles.dump'

    
foldername = ofolder + 'Cycles/' 
if os.path.isdir(foldername):
	shutil.rmtree(foldername,ignore_errors=True)

srcdata = ofolder + 'Cycle_10000/'
dstdata = ofolder + 'Cycles/Cycle_10000/'
shutil.copytree(srcdata,dstdata)

print('Dynamic isomerization started\n')

for Cycle_num in range(Ncycles1,Ncycles):
	print('Cycle Number = ',Cycle_num-Ncycles1+1, '   started')
	foldername = ofolder + 'Cycles/Cycle_' + str(Cycle_num) 
	if os.path.isdir(foldername):
		shutil.rmtree(foldername,ignore_errors=True)
	os.mkdir(foldername)
	os.chdir(foldername)
	
	dstdatar = ofolder + 'Cycles/Cycle_' + str(Cycle_num-1) + '/' + mdname + '_t' + str(1) + '.data'
	dstdataw = ofolder + 'Cycles/Cycle_' + str(Cycle_num) + '/' + mdname + '_t' + str(0) + '.data'
	
	datafile = open(dstdatar,'r')
	datawfile = open(dstdataw,'w')
	
	fname_dihed =  ofolder + 'Cycles/Cycle_'+ str(Cycle_num-1) + '/dump/dihed_' + str(dihedic)  + '.dat'
	Dihed_data_dummy = np.genfromtxt(fname_dihed,skip_header=9,skip_footer=0)
	
	Dihed_datateq = Dihed_data_dummy[Dihed_data_dummy[:,4]==81.0,:]
	Dihed_datatc = Dihed_data_dummy[Dihed_data_dummy[:,4]==108.0,:]
	Dihed_datact = Dihed_data_dummy[Dihed_data_dummy[:,4]==109.0,:]
	
	Dihed_data = np.vstack([Dihed_datateq,Dihed_datatc])
	Dihed_data = np.vstack([Dihed_data,Dihed_datact])
	
	count = 0
	count_trans_now = 0
	count_cis_now = 0
	
	count_trans_trans = 0
	count_cis_cis = 0
	
	count_trans_cis = 0
	count_cis_trans = 0
	for i in range(0,TNlines):
		count = count + 1
		line = datafile.readline()
		# print(line)
		if ((count > Ndihedline1 + 1)and(count<Ndihedline2)):
			# print(line)
			line1 = line.split(' ')
			# print(line1)
			if (int(line1[1])==81)or(int(line1[1])==108)or(int(line1[1])==109):
				dihedid = int(line1[2])
				# print(line1)
				dihedangle = 200
				for di in range(0,Dihed_data.shape[0]):
					if (int(Dihed_data[di,0])==dihedid):
						# print(Dihed_data[di,:])
						dihedangle = abs(Dihed_data[di,5])
				if (int(line1[1]) == 81): #existing equilibrium coefficients so change with prob
					if ((dihedangle>=90)and(dihedangle<=180)):
						count_trans_now = count_trans_now + 1
						if (np.random.uniform()<=probtc):
							line1[1] = str(108)
							count_trans_cis = count_trans_cis + 1
						else:
							line1[1] = str(109)
							count_trans_trans = count_trans_trans + 1
							
					elif ((dihedangle>=0)and(dihedangle<=90)):
						count_cis_now = count_cis_now + 1
						if (np.random.uniform()<=probct):
							line1[1] = str(109)
							count_cis_trans = count_cis_trans + 1
						else:
							line1[1] = str(108)	
							count_cis_cis = count_cis_cis + 1		
								
				elif (int(line1[1]) == 108): #current coefficients are trans - cis -- so if still in trans leave it
					if ((dihedangle>=90)and(dihedangle<=180)): #still in trans
						count_trans_now = count_trans_now + 1
						line1[1] = str(108)
						count_trans_cis = count_trans_cis + 1
					elif ((dihedangle>=0)and(dihedangle<=90)): #converted to cis
						count_cis_now = count_cis_now + 1
						if (np.random.uniform()<=probct):
							line1[1] = str(109)
							count_cis_trans = count_cis_trans + 1
						else:
							line1[1] = str(108)	
							count_cis_cis = count_cis_cis + 1	
							
				elif (int(line1[1]) == 109): #current coefficients are cis - trans -- so if still in cis leave it
					if ((dihedangle>=0)and(dihedangle<=90)): #still in cis
						count_cis_now = count_cis_now + 1
						line1[1] = str(109)
						count_cis_trans = count_cis_trans + 1
					elif ((dihedangle>=90)and(dihedangle<=180)): #converted to trans
						count_trans_now = count_trans_now + 1
						if (np.random.uniform()<=probtc):
							line1[1] = str(108)
							count_trans_cis = count_trans_cis + 1
						else:
							line1[1] = str(109)
							count_trans_trans = count_trans_trans + 1
							
			line = ''
			for j in range(0,5):
				line = line + line1[j] + ' '
			line = line + line1[5]
			# print(line) 
			datawfile.write(line)
		else:
			datawfile.write(line)
			
	datafile.close()
	datawfile.close()
	

	print('trans ', count_trans_now, ' cis ', count_cis_now)
	print('trans-cis ', count_trans_cis, ' trans-trans ', count_trans_trans)
	print('cis-trans ', count_cis_trans, ' cis-cis ', count_cis_cis)
		
	srcdata = ofolder + 'Cycles/Cycle_' + str(Cycle_num-1) + '/' + 'lammps2.in'
	dstdata = ofolder + 'Cycles/Cycle_' + str(Cycle_num) + '/' + 'lammps2.in'
	shutil.copy2(srcdata,dstdata)	
	
	srcdata = ofolder + 'Cycles/Cycle_' + str(Cycle_num-1) + '/' + 'my_hosts'
	dstdata = ofolder + 'Cycles/Cycle_' + str(Cycle_num) + '/' + 'my_hosts'
	shutil.copy2(srcdata,dstdata)	
	
	foldernamel = foldername + '/log' 
	os.mkdir(foldernamel)

	foldernamed = foldername + '/dump' 
	os.mkdir(foldernamed)
	
	subprocess.call(lammpscmd, shell=True)
	os.chdir(ofolder)
	pristr = 'Cycle Number = %7d out of %7d = %4.2f %% done :)' % (Cycle_num-Ncycles1+1, Dy_Ncycles, (Cycle_num-Ncycles1+1)*100.0/Dy_Ncycles)
	print(pristr)
	print('Time taken till now ', datetime.now()-starttime,'\n')

print('All done :)')
