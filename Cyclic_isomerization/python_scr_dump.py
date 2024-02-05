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
import os
import shutil

ofolder= os.getcwd()
ofolder = ofolder + '/'
mdname = 'M10D1A1_C25'

dstdataw = ofolder + 'All_cycles.dump'
datawfile = open(dstdataw,'w')


NNN = 5
Ndump = 10000
Nlineseach = 7287
TNlines = Nlineseach*201
foldername = ofolder + 'Cycles/Cycle_1' 
	
dstdatar = foldername + '/excite.dump'
	
datafile = open(dstdatar,'r')
line1 = datafile.readline()
datafile.close()


count = 0
for Cycle_num in range(1,2):
	foldername = ofolder + 'Cycles/Cycle_' + str(Cycle_num) 
	
	dstdatar = foldername + '/excite.dump'
	
	datafile = open(dstdatar,'r')
	

	
	
	for i in range(0,TNlines):
		
		line = datafile.readline()
		if (i>=0):
			if (line==line1):
				datawfile.write(line)
				datawfile.write(str(count*Ndump) + line1[-1])
				line = datafile.readline()
				count = count + 1
			else:
				datawfile.write(line)
				

	datafile.close()
	
	os.chdir(ofolder)
	print('Cycle Number =',Cycle_num)
	

for Cycle_num in range(2,NNN):
	foldername = ofolder + 'Cycles/Cycle_' + str(Cycle_num) 
	
	dstdatar = foldername + '/excite.dump'
	
	datafile = open(dstdatar,'r')
	

	
	
	for i in range(0,TNlines):
		
		line = datafile.readline()
		if (i>=Nlineseach):
			if (line==line1):
				datawfile.write(line)
				datawfile.write(str(count*Ndump) + line1[-1])
				line = datafile.readline()
				count = count + 1
			else:
				datawfile.write(line)
				

	datafile.close()
	
	os.chdir(ofolder)
	print('Cycle Number =',Cycle_num)
	
datawfile.close()
