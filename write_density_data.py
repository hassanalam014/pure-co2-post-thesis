# Date: April 2017
#
# Description: The purpose of this file is to write experimental data files, in CSV format,
# 			   based on the Helmholtz free energy proposed in the paper:
#			   Wu et al, 2011. Journal of Physical and Chemical Reference Data (40).
#			   #Valid from the triple point (2.2Pa, 131.66K) to 550K and up to 50MPa.
#			   Not valid near the critical point (5.3368MPa, 400.378K).
#

import os,sys,math,matplotlib.pyplot as plt,numpy as npy
import csv
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
# from DME_EOS import *
from CO2_EOS import *

#================================================================
# PARAMETERS
#================================================================

#Setting whether or not to plot results.
plot_result = True

#Setting the name of the output folder.
output_folder = 'Data'
#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)

#Initializing the array of densities.
R0 = npy.linspace(0.03,0.35,50)

#Molar mass of Dimethyl Ether.
M = 46.07

# Array of temperatures to create isotherms.
Tlist = [423,433,443,453,463,473,483]

#================================================================
# CALCULATING THEORETICAL ISOTHERMS AND WRITING OUTPUT TO FILE
#================================================================

length = len(R0)+6

for i in range(0,len(Tlist)):
	
	exec "T%s_P = calculatePressureDME(float(%s),R0)" % (Tlist[i],Tlist[i])
	exec "P0 = T%s_P" % (Tlist[i])
	
	output_array = npy.zeros((length,4))
	
	output_array[4][0] = len(R0)
	
	for j in range(0,len(R0)):
		output_array[j+6][0] = P0[j]
		output_array[j+6][1] = float(Tlist[i])
		output_array[j+6][2] = R0[j]
		output_array[j+6][3] = M
	
	exec "name = './Data/%sK_DME_PVT.csv'" % (Tlist[i])
	with open(name, 'wb') as f:
		writer = csv.writer(f)
		writer.writerows(output_array)

#================================================================
# PLOTTING RESULTS
#================================================================

if plot_result:
	plt.figure(num=None, figsize=(12, 10), dpi=80, facecolor='w', edgecolor='k')
	ax = plt.axes()
	plt.plot(T423_P,R0,'-r')
	plt.plot(T463_P,R0,'-g')
	plt.plot(T483_P,R0,'-b')

	# plt.plot(T423_P,R0,'or')
	# plt.plot(T463_P,R0,'sg')
	# plt.plot(T483_P,R0,'^b')

	plt.xlabel('Pressure P (MPa)',fontsize=14)
	plt.ylabel(r'Density $\rho$ ($g/cm^3$)',fontsize=14)
	plt.axis([0,85,0.0,0.38])

	plt.show()
