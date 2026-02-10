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

from matplotlib.ticker import AutoMinorLocator
from loadExperimentalData import *
from findVectors import findVectors

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
R = npy.linspace(0.03,0.35,10)

#Molar mass of Dimethyl Ether.
M = 46.07

# Array of temperatures to create isotherms.
Tlist = [323,373,380,400,420,490,660]

#================================================================
# CALCULATING THEORETICAL ISOTHERMS AND WRITING OUTPUT TO FILE
#================================================================

# length = len(R0)+6

# for i in range(0,len(Tlist)):
	
	# exec "T%s_P = calculatePressureDME(float(%s),R0)" % (Tlist[i],Tlist[i])
	# exec "P0 = T%s_P" % (Tlist[i])
	
	# exec "T%s_P = result[0]" % (Tlist[i])

	# exec "vector_%s = findVectors(T%s_P,R0,P0_%s,R0_%s)" % (Tlist[i],Tlist[i],Tlist[i],Tlist[i])

	# output_array = npy.zeros((length,4))
	
	# output_array[4][0] = len(R0)
	
	# for j in range(0,len(R0)):
	# 	output_array[j+6][0] = P0[j]
	# 	output_array[j+6][1] = float(Tlist[i])
	# 	output_array[j+6][2] = R0[j]
	# 	output_array[j+6][3] = M
	
	# exec "name = './Data/%sK_DME_PVT.csv'" % (Tlist[i])
	# with open(name, 'wb') as f:
	# 	writer = csv.writer(f)
	# 	writer.writerows(output_array)

#================================================================
# PLOTTING RESULTS
#================================================================

arrow_ls = 'dashdot'
show_arrows = True

#Setting font size
axis_size = 20
title_size = 20
size = 14
label_size = 20
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

#Setting saved image properties
img_extension = '.pdf'
img_dpi = None
output_folder = 'plot_density'

#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)

#==================================================================================
#P versus R plots.
figPVSR = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()
'''
#plt.plot(P0_311,R0_311,'ob')
plt.plot(P0_323,R0_323,'og')
#plt.plot(P0_340A,R0_340A,'ok')
plt.plot(P0_373,R0_373,'sr')
plt.plot(P0_380,R0_380,'^c')
#plt.plot(P0_380A,R0_380A,'oc')
#plt.plot(P0_381,R0_381,'ok')
#plt.plot(P0_388,R0_388,'oc')
plt.plot(P0_400,R0_400,'vm')
#plt.plot(P0_406,R0_406,'om')
#plt.plot(P0_410,R0_410,'ok')
plt.plot(P0_420,R0_420,'*y')
#plt.plot(P0_430,R0_430,'og')
plt.plot(P0_490,R0_490,'dk')
#plt.plot(P0_520,R0_520,'ob')
plt.plot(P0_660,R0_660,'>',color='orange')
# plt.plot(P0_1100,R0_1100,'<',color='gray')
#plt.plot(P0_425,R0_425,'oy')
#plt.plot(P0_307,R0_307,'dk')
#plt.plot(P0_310,R0_310,'dy')
#plt.plot(P0_313,R0_313,'dk')
#plt.plot(P0_320,R0_320,'dk')
#plt.plot(T311_P,x0,'b')
# plt.plot(T323_P,x0,'g')
# plt.plot(T373_P,x0,'r')
#plt.plot(T340_P,x0,'k')
# plt.plot(T380_P,x0,'c')
#plt.plot(T381_P,x0,'k')
#plt.plot(T388_P,x0,'c')
# plt.plot(T400_P,x0,'m')
#plt.plot(T406_P,x0,'m')
#plt.plot(T410_P,x0,'k')
# plt.plot(T420_P,x0,'y')
#plt.plot(T430_P,x0,'g')
# plt.plot(T490_P,x0,'k')
#plt.plot(T520_P,x0,'b')
# plt.plot(T660_P,x0,color='orange')
# plt.plot(T1100_P,x0,color='gray')
#plt.plot(T425_P,x0,'y')
#plt.plot(T307_P,x0)
#plt.plot(T310_P,x0,'y')
#plt.plot(T313_P,x0)
#plt.plot(T320_P,x0,'k')
'''
# if show_arrows:
	# for i in range(10,len(P0_323)):
		#ax.arrow(P0_323[i],R0_323[i],vector_323[0][i],vector_323[1][i],width=0.002,head_width=0.03,head_length=0.80,length_includes_head=True,fc='g',ec='g')
		#head_angle = math.atan2(vector_323[1][i],vector_323[0][i])
		#ax.arrow(P0_323[i],R0_323[i],vector_323[0][i],vector_323[1][i],head_length=abs(0.65*math.cos(head_angle)),head_width=abs(0.01*math.sin(head_angle)),length_includes_head=True,fc='g',ec='g')
		# ax.arrow(P0_323[i],R0_323[i],vector_323[0][i],vector_323[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='g',ec='g')
		#ax.annotate("",xy=(1.5, 1.5),xytext=(0, 0),arrowprops=dict(arrowstyle="->"))
	# for i in range(9,len(P0_373)):
	# 	ax.arrow(P0_373[i],R0_373[i],vector_373[0][i],vector_373[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='r',ec='r')
	# for i in range(0,len(P0_380)):
	# 	ax.arrow(P0_380[i],R0_380[i],vector_380[0][i],vector_380[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='c',ec='c')
	# for i in range(0,len(P0_400)):
	# 	ax.arrow(P0_400[i],R0_400[i],vector_400[0][i],vector_400[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='m',ec='m')
	# for i in range(0,len(P0_420)):
	# 	ax.arrow(P0_420[i],R0_420[i],vector_420[0][i],vector_420[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='y',ec='y')
	# for i in range(0,len(P0_490)):
	# 	ax.arrow(P0_490[i],R0_490[i],vector_490[0][i],vector_490[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='k',ec='k')
	# for i in range(0,len(P0_660)):
	# 	ax.arrow(P0_660[i],R0_660[i],vector_660[0][i],vector_660[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='orange',ec='orange')
	# for i in range(0,len(P0_1100)):
	# 	ax.arrow(P0_1100[i],R0_1100[i],vector_1100[0][i],vector_1100[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='grey',ec='grey')
plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel(r'Density $\rho$ ($g/cm^3$)',fontsize=axis_size)
plt.legend(('Garg 323K','Garg 373K','Xiong 380K','Xiong 400K','Xiong 420K','Angus 490K','Angus 660K','Angus 1100K'),loc=2,fontsize=size,numpoints=1)
#plt.title(title,fontsize=title_size)
plt.axis([0,68,0,1.0])
minorLocator = AutoMinorLocator()
ax.xaxis.set_minor_locator(minorLocator)
plt.tick_params(which='both', width=1)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4)
minorLocator = AutoMinorLocator()
ax.yaxis.set_minor_locator(minorLocator)
plt.tick_params(which='both', width=1)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4)
# figPVSR.savefig('./'+output_folder+r'\pure_CO2_density_isotherm_'+param_set+img_extension,dpi=img_dpi)
'''
# plt.figure(num=None, figsize=(12, 10), dpi=80, facecolor='w', edgecolor='k')
# ax = plt.axes()
plt.plot(T323_P,R0,'-r')
plt.plot(T373_P,R0,'-g')
plt.plot(T380_P,R0,'-b')
plt.plot(T400_P,R0,'-r')
plt.plot(T420_P,R0,'-g')
plt.plot(T490_P,R0,'-b')
plt.plot(T660_P,R0,'-r')
'''
plt.plot(P0,R0,'or')
# plt.plot(T463_P,R0,'sg')
# plt.plot(T483_P,R0,'^b')

plt.xlabel('Pressure P (MPa)',fontsize=14)
plt.ylabel(r'Density $\rho$ ($g/cm^3$)',fontsize=14)
# plt.axis([0,85,0.0,0.38])

plt.figure(num=None, figsize=(12, 10), dpi=80, facecolor='w', edgecolor='k')
ax = plt.axes()
plt.plot(T0,R0,'or')

plt.xlabel('Temperature T (K)',fontsize=14)
plt.ylabel(r'Density $\rho$ ($g/cm^3$)',fontsize=14)
# plt.axis([0,85,0.0,0.38])

plt.show()
