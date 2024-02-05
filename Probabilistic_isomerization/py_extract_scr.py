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

ofolder= os.getcwd()
ofolder = ofolder + '/'
mdname = 'M10D1A1_C25_eq2'
lammpscmd = '/usr/lib64/openmpi/bin/mpirun -np 8 --hostfile my_hosts /home/Akhil_CRSM/Installations/lammps/src/lmp_mpi -in lammps2.in >> log/terminal.log'
#subprocess.call(lammpscmd, shell=True)

#Dihedrals 81 is thermal, 108 - is tc and 109 is ct
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

srcdata = ofolder + 'Cycle_' + str(10000) + '/' + file1name
dstdata = ofolder + file1name
shutil.copy2(srcdata,dstdata)	
	
srcdata = ofolder + 'Cycle_' + str(10000) + '/' + file2name
dstdata = ofolder + file2name
shutil.copy2(srcdata,dstdata)	

srcdata = ofolder + 'Cycle_' + str(10000) + '/' + dumpfilea
dstdata = ofolder + dumpfilea
shutil.copy2(srcdata,dstdata)	

print('Extraction for dynamic isomerization extraction started\n')

retry = 5
sleeptime = 180.0

A0=np.genfromtxt(file1name)
LastIter = A0[-1,0]
LastTime = A0[-1,1]

time.sleep(0.0)
for Cycle_num in range(Ncycles1,Ncycles):
	
	filetocheck =  ofolder + 'Cycles/Cycle_' + str(Cycle_num) + '/' + mdname + '_t' + str(1) + '.data'
	retries = 0
	fileexist = 0 
	print('Extracting for Cycle Number = %7d started' % (Cycle_num - Ncycles1 + 1))		
	while ((retries < retry)and(fileexist==0)):
		if(os.path.exists(filetocheck)):
			fileexist = 1
			time.sleep(0.5)
			print('File exists started extracting')
			filename = 'Cycles/Cycle_' + str(Cycle_num) +'/log/l_npt_ex.log'
			A1=np.genfromtxt(filename,skip_header=27,skip_footer=38)

			dumpname = 'Cycles/Cycle_' + str(Cycle_num) + '/excite.dump'

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

			Dir = np.zeros((Niters,4))
			while (count < Nlines):
				count += 1
				line = dumpfile.readline()
				if (line == line1):
					count_iter +=1
					line_tstep = dumpfile.readline()
					Tstepc = int(line_tstep)
					Dir[count_iter-1,0] = Tstepc
					# print(line,count,count_iter,Tstep)	
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
						# print(line.split(' '))
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

			# print(A[0,0,:])
			for i in range(1,Niters):	
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
				
				#print(Ndir)
				C = np.zeros((Ndir,21))
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
					

					C[j,12] = C[j,8]/C[j,11]
					C[j,13] = C[j,9]/C[j,11]
					C[j,14] = C[j,10]/C[j,11]
					
					C[j,15] = (3.0*C[j,12]**2-1.0)/2.0
					C[j,16] = (3.0*C[j,13]**2-1.0)/2.0
					C[j,17] = (3.0*C[j,14]**2-1.0)/2.0	

				Dir[i,1] = np.mean(C[:,15])
				Dir[i,2] = np.mean(C[:,16])
				Dir[i,3] = np.mean(C[:,17])

			Dir=Dir[Dir[:,1]>0,:]
			A1=np.hstack([A1,Dir[:,1:4]])
			A1[:,0] = A1[:,0] + LastIter
			A1[:,1] = A1[:,1] + LastTime
			LastIter = A1[-1,0]
			LastTime = A1[-1,1]
			file1name = 'Data_TTPDLO.dat'
			f1=open(file1name,'a')
			np.savetxt(f1,A1,fmt='%15.6f')
			f1.close()

			B1 = np.zeros([Niters-1,Nisomers + 6])
			icount = 0

			for i in range(1,Niters):
				
				
				fname =  'Cycles/Cycle_'+ str(Cycle_num) + '/dump/dihed_' + str(i*NDiheddump) + '.dat'
				C1=np.genfromtxt(fname,skip_header=9,skip_footer=0)
				
				C2a=C1[C1[:,4]==81.0,:]
				C2b=C1[C1[:,4]==108.0,:]
				C2c=C1[C1[:,4]==109.0,:]


				if ((len(C2a)>0)and(len(C2b)>0)):
					C2d = np.vstack((C2a,C2b))
					if (len(C2c>0)):
						C2=np.vstack((C2d,C2c))
					else:
						C2 = C2d
				elif ((len(C2a)>0)and(len(C2b)==0)):
					C2d = C2a
					if (len(C2c>0)):
						C2=np.vstack((C2d,C2c))
					else:
						C2 = C2d			
				elif ((len(C2a)==0)and(len(C2b)>0)):
					C2d = C2b
					if (len(C2c>0)):
						C2=np.vstack((C2d,C2c))
					else:
						C2 = C2d			
				elif ((len(C2a)==0)and(len(C2b)==0)):
					if (len(C2c>0)):
						C2=C2c
					else:
						C2 = []	
				C2 = C2[C2[:,0].argsort()]	
				#C2 = np.vstack([C2a,C2b])
				for j in range(0,Nisomers):
					B1[i-1,j+2] = C2[j,5]				
			
				C2[:,5] = abs(C2[:,5])
				D1 = np.zeros([C2.shape[0],3])
				# print(C2[:,0:6])
				
				# print(C2[:,0:6])
				if (Actual_iter>0):
					if (np.sum(C2[:,5]>=90.0) > 0):
						D1[C2[:,5]>=90.0,0] = 1.0
					if (np.sum(C2[:,5]<90.0) > 0):
						D1[C2[:,5]<90.0 ,0] = 0.0
					if (np.sum(Cold[:,5]>=90.0) > 0):
						D1[Cold[:,5]>=90.0,1] = 1.0
					if (np.sum(Cold[:,5]<90.0) > 0):
						D1[Cold[:,5]<90.0,1] = 0.0
					
					D1[:,2] = np.absolute(D1[:,0] - D1[:,1])

				Cold = C2
				B1[i-1,0] = A1[i-1,0]
				B1[i-1,1] = A1[i-1,1]

				B1[i-1,Nisomers + 2] = np.mean(C2[:,5])
				B1[i-1,Nisomers + 3] = np.sum(C2[:,5]>=90.0)
				B1[i-1,Nisomers + 4] = C2.shape[0] - B1[i-1,Nisomers + 3]
				B1[i-1,Nisomers + 5] = np.sum(D1[:,2])
				icount = icount + 1
				Actual_iter =Actual_iter + 1

			file2name = 'Data_Diheds.dat'
			f2=open(file2name,'a')
			np.savetxt(f2,B1,fmt='%15.6f')
			f2.close()
			
			pristr = 'Cycle Number = %7d out of %7d = %4.2f %% done :)\n' % (Cycle_num-Ncycles1+1, Dy_Ncycles, (Cycle_num-Ncycles1+1)*100.0/Dy_Ncycles)
			print(pristr)
			
			if ((Cycle_num % Nmesogendump) ==0):
				datawfile = open(dumpfilea,'a')

				foldername = ofolder + 'Cycles/Cycle_' + str(Cycle_num) 
				dstdatar = foldername + '/excite.dump'
				datafile = open(dstdatar,'r')
				line1 = datafile.readline()
				for i in range(0,TNlinesdump):
					line = datafile.readline()
					if (i>=Nlineseach*dihedic-1):
						# print(line)
						if (line==line1):
							datawfile.write(line)
							datawfile.write(str(Cycle_num) + line1[-1])
							line = datafile.readline()
							count = count + 1
						else:
							datawfile.write(line)
				datafile.close()
				datawfile.close()
				
				print('Updated the dump file with cycle ', Cycle_num ,'\n')	
			
			if ((Cycle_num % Ncycledelete)==0):
				for cnum in range(Cycle_num - 2*Ncycledelete +1 ,Cycle_num - Ncycledelete + 1):
					oldfolder = ofolder + 'Cycles/Cycle_' + str(cnum) + '/' 
					# print(oldfolder)
					if os.path.isdir(oldfolder):
						shutil.rmtree(oldfolder,ignore_errors=True)
				print('Cycles',Cycle_num - 2*Ncycledelete + 1,' to ',Cycle_num - Ncycledelete, ' are deleted now\n')				
			
		else:	
			fileexist = 0
			print('File doesnt exist so sleeping out for ', sleeptime, ' seconds')
			time.sleep(sleeptime)
			retries = retries + 1

print('All done :)')
