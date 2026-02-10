# Date: April 2017
#
# Description: The purpose of this file is to plot Carbon Dioxide (CO2) density information based on experiment and theory for comparison.
#

import os,sys,math,matplotlib.pyplot as plt,numpy as npy
from matplotlib.ticker import AutoMinorLocator
from all_s_params import *
from loadExperimentalData import *
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
from loadPhysicalConstants import *
from P4D_plotting_parameters import *
from findVectors import findVectors
from wrapperFunctions import calculatePressure,calculateTemperature,calculateDensity,calculateIsothermalCompressibility,calculateThermalExpansivity,calculateIsobaricHeatCapacity,calculateIsochoricHeatCapacity,calculateEntropy,calculateSecondVirialCoefficient,calculatePureCoexistence
from calculatePureVariables import calculateNewMolecularParameters,calculateCharacteristicParametersGamma,calculatePureCriticalPoint

#CO2 Molecular weight
M = 44.01

#Setting which set of parameters to use for calculation.
param_set = 'Self'

#Initializing the array of densities.
x0 = npy.linspace(0.01,1.2,300)
y0 = npy.linspace(0.01,1.2,300)

if param_set == 'Kilpatrick':
	title = 'Density of CO2 using Kilpatrick 1986 SL Parameters'
	Pstar = Kilpatrick_Pstar
	Tstar = Kilpatrick_Tstar
	Rstar = Kilpatrick_Rstar
elif param_set == 'Kiszka':
	title = 'Density of CO2 using Kiszka 1988 SL Parameters'
	Pstar = Kiszka_Pstar
	Tstar = Kiszka_Tstar
	Rstar = Kiszka_Rstar
elif param_set == 'Pope':
	title = 'Density of CO2 using Pope 1991 SL Parameters'
	Pstar = Pope_Pstar
	Tstar = Pope_Tstar
	Rstar = Pope_Rstar
elif param_set == 'Hariharan':
	title = 'Density of CO2 using Hariharan 1993 SL Parameters'
	Pstar = Hariharan_Pstar
	Tstar = Hariharan_Tstar
	Rstar = Hariharan_Rstar
elif param_set == 'Garg':
	title = 'Density of CO2 using Garg 1994 SL Parameters'
	Pstar = Garg_Pstar
	Tstar = Garg_Tstar
	Rstar = Garg_Rstar
elif param_set == 'Xiong':
	title = 'Density of CO2 using Xiong 1995 SL Parameters'
	Pstar = Xiong_Pstar
	Tstar = Xiong_Tstar
	Rstar = Xiong_Rstar
elif param_set == 'Cao':
	title = 'Density of CO2 using Cao 2010 SL Parameters'
	Pstar = Cao_Pstar
	Tstar = Cao_Tstar
	Rstar = Cao_Rstar
elif param_set == 'Nalawade':
	title = 'Density of CO2 using Nalawade 2006 SL Parameters'
	Pstar = Nalawade_Pstar
	Tstar = Nalawade_Tstar
	Rstar = Nalawade_Rstar
elif param_set == 'Self':
	title = 'Density of CO2 using our own SL Parameters'
	Pstar = Self_Pstar
	Tstar = Self_Tstar
	Rstar = Self_Rstar
elif param_set == 'Funami':
	title = 'Density of CO2 using Funami 2007 SL Parameters'
	Pstar = Funami_Pstar
	Tstar = Funami_Tstar
	Rstar = Funami_Rstar
elif param_set == 'Arce':
	title = 'Density of CO2 using Arce 2009 SL Parameters'
	Pstar = Arce_Pstar
	Tstar = Arce_Tstar
	Rstar = Arce_Rstar
elif param_set == 'Doghieri':
	title = 'Density of CO2 using Doghieri 1996 SL Parameters'
	Pstar = Doghieri_Pstar
	Tstar = Doghieri_Tstar
	Rstar = Doghieri_Rstar

gamma, vh, epsilon = calculateNewMolecularParameters(Pstar,Tstar,Rstar,M)
vh = vh/NA
epsilon = epsilon/NA
print('gamma = {}, vh = {}, and epsilon = {}'.format(gamma,vh,epsilon))

if param_set == 'Wang':
	title = 'Density of CO2 using Wang 1991 SL Parameters'
	
	result = calculatePressure(311.85,Wang_Pstar,278.518,Wang_Rstar,M,x0)
	T311_P = result[0]

	result = calculatePressure(323.15,Wang_Pstar,278.280,Wang_Rstar,M,x0)
	T323_P = result[0]

	result = calculatePressure(373.15,Wang_Pstar,274.150,Wang_Rstar,M,x0)
	T373_P = result[0]

	result = calculatePressure(380.0,Wang_Pstar,274.154,Wang_Rstar,M,x0)
	T380_P = result[0]

	result = calculatePressure(400.0,Wang_Pstar,271.540,Wang_Rstar,M,x0)
	T400_P = result[0]

	result = calculatePressure(420.0,Wang_Pstar,286.322,Wang_Rstar,M,x0)
	T420_P = result[0]
	
	result = calculatePressure(280.0,Wang_Pstar,278.15,Wang_Rstar,M,x0)
	T280_P = result[0]
else:
	#==============================================================================================================
	#Calculating Isotherms.
	#==============================================================================================================
	
	result = calculatePressure(311.85,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T311_P = result[0]

	result = calculatePressure(323.15,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T323_P = result[0]
	vector_323 = findVectors(T323_P,x0,P0_323,R0_323)

	result = calculatePressure(373.15,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T373_P = result[0]
	vector_373 = findVectors(T373_P,x0,P0_373,R0_373)
	result = calculateIsothermalCompressibility(T373_P,373.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T373_B = result[1]

	result = calculatePressure(380.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T380_P = result[0]
	vector_380 = findVectors(T380_P,x0,P0_380,R0_380)
	result = calculateIsothermalCompressibility(T380_P,380.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T380_B = result[1]

	result = calculatePressure(400.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T400_P = result[0]
	vector_400 = findVectors(T400_P,x0,P0_400,R0_400)
	result = calculateIsothermalCompressibility(T400_P,400.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T400_B = result[1]
	
	result = calculatePressure(410.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T410_P = result[0]

	result = calculatePressure(420.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T420_P = result[0]
	vector_420 = findVectors(T420_P,x0,P0_420,R0_420)
	
	result = calculatePressure(430.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T430_P = result[0]
	result = calculateIsothermalCompressibility(T430_P,430.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T430_B = result[1]
	
	result = calculatePressure(490.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T490_P = result[0]
	vector_490 = findVectors(T490_P,x0,P0_490,R0_490)
	result = calculateIsothermalCompressibility(T490_P,490.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T490_B = result[1]
	
	result = calculatePressure(520.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T520_P = result[0]
	result = calculateIsothermalCompressibility(T520_P,520.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T520_B = result[1]
	
	result = calculatePressure(660.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T660_P = result[0]
	vector_660= findVectors(T660_P,x0,P0_660,R0_660)
	result = calculateIsothermalCompressibility(T660_P,660.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T660_B = result[1]
	
	result = calculatePressure(1100.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T1100_P = result[0]
	vector_1100= findVectors(T1100_P,x0,P0_1100,R0_1100)
	result = calculateIsothermalCompressibility(T1100_P,1100.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T1100_B = result[1]
	
	result = calculatePressure(381.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T381_P = result[0]
	
	result = calculatePressure(388.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T388_P = result[0]
	
	result = calculatePressure(406.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T406_P = result[0]
	
	result = calculatePressure(425.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T425_P = result[0]
	
	result = calculatePressure(307.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T307_P = result[0]
	
	result = calculatePressure(310.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T310_P = result[0]
	
	result = calculatePressure(313.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T313_P = result[0]
	
	result = calculatePressure(320.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T320_P = result[0]
	
	result = calculatePressure(340.0,x0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	T340_P = result[0]
	
	#==============================================================================================================
	
	T_SVC = npy.linspace(100,1100,100)
	result = calculateSecondVirialCoefficient(T_SVC,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P0MPa_BC = result[1]
	
	#==============================================================================================================
	# Calculating Isobars.
	#==============================================================================================================
	
	result = calculateTemperature(0.5,y0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P0pt5MPa_T = result[0]
	
	result = calculateTemperature(1.0,y0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P1MPa_T = result[0]
	
	result = calculateTemperature(6.0,y0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P6MPa_T = result[0]
	
	result = calculateTemperature(12.0,y0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P12MPa_T = result[0]
	
	result = calculateTemperature(16.0,y0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P16MPa_T = result[0]
	result = calculateThermalExpansivity(16.0,P16MPa_T,y0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P16MPa_A = result[1]
	
	result = calculateTemperature(30.0,y0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P30MPa_T = result[0]
	vector_30MPa = findVectors(P30MPa_T,y0,T0_30MPa,R0_30MPa)
	result = calculateThermalExpansivity(30.0,P30MPa_T,y0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P30MPa_A = result[1]
	
	result = calculateTemperature(50.0,y0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P50MPa_T = result[0]
	vector_50MPa = findVectors(P50MPa_T,y0,T0_50MPa,R0_50MPa)
	result = calculateThermalExpansivity(50.0,P50MPa_T,y0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P50MPa_A = result[1]
	result = calculateIsobaricHeatCapacity(50.0,P50MPa_T,y0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P50MPa_CP = result[1]
	result = calculateEntropy(50.0,P50MPa_T,y0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P50MPa_S = result[1]
	
	result = calculateTemperature(65.0,y0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P65MPa_T = result[0]
	vector_65MPa = findVectors(P65MPa_T,y0,T0_65MPa,R0_65MPa)
	result = calculateThermalExpansivity(65.0,P65MPa_T,y0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P65MPa_A = result[1]
	result = calculateIsobaricHeatCapacity(65.0,P65MPa_T,y0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P65MPa_CP = result[1]
	result = calculateEntropy(65.0,P65MPa_T,y0,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P65MPa_S = result[1]
	
	result = calculateIsochoricHeatCapacity(P0_CV1,T0_CV1,0.0550125,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	R55_CV = result[1]
	
	result = calculateIsochoricHeatCapacity(P0_CV2,T0_CV2,1.0232325,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	R1023_CV = result[1]
	
	result = calculateIsochoricHeatCapacity(P0_CV3,T0_CV3,0.79218,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	R792_CV = result[1]
	
	#==============================================================================================================
	#Calculating Pressure Deviations.
	#==============================================================================================================
	
	result = calculatePureCoexistence(T0_Coex_Duschek,M0_Coex_Duschek,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
	P_dev_D = result[1]
	R_dev_G_D = result[2]
	R_dev_L_D = result[3]
	deviation_L_D = range(0,len(P_dev_D))
	deviation_G_D = range(0,len(P_dev_D))
	for i in range(0,len(P_dev_D)):
		deviation_L_D[i] = 100*abs(R_dev_L_D[i]-R0_CoeL_Duschek[i])/R0_CoeL_Duschek[i]
		deviation_G_D[i] = 100*abs(R_dev_G_D[i]-R0_CoeG_Duschek[i])/R0_CoeG_Duschek[i]
		#print deviation_G_D[i]
		#print R_dev_G_D[i]
		#print R0_CoeG_Duschek[i]
	
	temps = [323,373,380,400,420,490,660,1100]
	for k in range(0,len(temps)):
		exec "result = calculateDensity(P0_%s,T0_%s,M0_%s,Pstar=Kilpatrick_Pstar,Tstar=Kilpatrick_Tstar,Rstar=Kilpatrick_Rstar)" % (temps[k],temps[k],temps[k])
		exec "R_dev_%s = result[1]" % (temps[k])
		exec "Kilpatrick_deviation_%s = range(0,len(P0_%s))" % (temps[k],temps[k])
		exec "for i in range(0,len(P0_%s)):\n	Kilpatrick_deviation_%s[i] = 100*abs(R_dev_%s[i]-R0_%s[i])/R0_%s[i]" % (temps[k],temps[k],temps[k],temps[k],temps[k])
	for k in range(0,len(temps)):
		exec "result = calculateDensity(P0_%s,T0_%s,M0_%s,Pstar=Self_Pstar,Tstar=Self_Tstar,Rstar=Self_Rstar)" % (temps[k],temps[k],temps[k])
		exec "R_dev_%s = result[1]" % (temps[k])
		exec "Self_deviation_%s = range(0,len(P0_%s))" % (temps[k],temps[k])
		exec "for i in range(0,len(P0_%s)):\n	Self_deviation_%s[i] = 100*abs(R_dev_%s[i]-R0_%s[i])/R0_%s[i]" % (temps[k],temps[k],temps[k],temps[k],temps[k])
	
	press = ['6MPa','12MPa','16MPa','30MPa','50MPa','65MPa']
	for k in range(0,len(press)):
		exec "result = calculateDensity(P0_%s,T0_%s,M0_%s,Pstar=Kilpatrick_Pstar,Tstar=Kilpatrick_Tstar,Rstar=Kilpatrick_Rstar)" % (press[k],press[k],press[k])
		exec "R_dev_%s = result[1]" % (press[k])
		exec "Kilpatrick_deviation_%s = range(0,len(P0_%s))" % (press[k],press[k])
		exec "for i in range(0,len(P0_%s)):\n	Kilpatrick_deviation_%s[i] = 100*abs(R_dev_%s[i]-R0_%s[i])/R0_%s[i]" % (press[k],press[k],press[k],press[k],press[k])
	for k in range(0,len(press)):
		exec "result = calculateDensity(P0_%s,T0_%s,M0_%s,Pstar=Self_Pstar,Tstar=Self_Tstar,Rstar=Self_Rstar)" % (press[k],press[k],press[k])
		exec "R_dev_%s = result[1]" % (press[k])
		exec "Self_deviation_%s = range(0,len(P0_%s))" % (press[k],press[k])
		exec "for i in range(0,len(P0_%s)):\n	Self_deviation_%s[i] = 100*abs(R_dev_%s[i]-R0_%s[i])/R0_%s[i]" % (press[k],press[k],press[k],press[k],press[k])
	
	group = ['Kilpatrick','Kiszka','Wang','Pope','Hariharan','Garg','Xiong','Cao','Nalawade','Self','Funami','Arce','Doghieri']
	temp = 490
	for k in range(0,len(group)):
		exec "result = calculateDensity(P0_%s,T0_%s,M0_%s,Pstar=%s_Pstar,Tstar=%s_Tstar,Rstar=%s_Rstar)" % (temp,temp,temp,group[k],group[k],group[k])
		exec "R_dev_%s = result[1]" % (group[k])
		exec "%s_deviation = range(0,len(P0_%s))" % (group[k],temp)
		exec "for i in range(0,len(P0_%s)):\n	%s_deviation[i] = 100*abs(R_dev_%s[i]-R0_%s[i])/R0_%s[i]" % (temp,group[k],group[k],temp,temp)
	
	#==============================================================================================================

if param_set != 'Wang':
	r = (Pstar*M)/(8.31*Tstar*Rstar)
	Pc = Pstar*(2/(1+math.sqrt(r))**2)*(r*math.log(1.0+1.0/math.sqrt(r))+0.5-math.sqrt(r))
	Tc = Tstar*(2*r/(1+math.sqrt(r))**2)
	print('The SL method critical temperature predicted by {} is {} K at {} MPa.'.format(param_set,Tc,Pc))

#if param_set != 'Wang':
	#result = calculatePureCriticalPoint(Pstar,Tstar,Rstar,Ms)
	#Pc = result[0]
	#Tc = result[1]
	#print('The calculated critical temperature predicted by {} is {} K at {} MPa.'.format(param_set,Tc,Pc))

T0_Coex = npy.concatenate((T0_Coex_Newitt,T0_Coex_Angus,T0_Coex_Duschek),axis=0)
T0_Theory = npy.linspace(min(T0_Coex),0.99*max(T0_Coex),100)
result = calculatePureCoexistence(T0_Theory,M,Pstar=Pstar,Tstar=Tstar,Rstar=Rstar)
P0_Theory = result[1]
rho_g = result[2]
rho_l = result[3]
H_vap = result[4]

arrow_ls = 'dashdot'
show_arrows = True

#==================================================================================
#P versus R plots.
figPVSR = plt.figure(num=None, figsize=(img_vert, img_hori), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()
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
plt.plot(P0_1100,R0_1100,'<',color='gray')
#plt.plot(P0_425,R0_425,'oy')
#plt.plot(P0_307,R0_307,'dk')
#plt.plot(P0_310,R0_310,'dy')
#plt.plot(P0_313,R0_313,'dk')
#plt.plot(P0_320,R0_320,'dk')
#plt.plot(T311_P,x0,'b')
plt.plot(T323_P,x0,'g')
plt.plot(T373_P,x0,'r')
#plt.plot(T340_P,x0,'k')
plt.plot(T380_P,x0,'c')
#plt.plot(T381_P,x0,'k')
#plt.plot(T388_P,x0,'c')
plt.plot(T400_P,x0,'m')
#plt.plot(T406_P,x0,'m')
#plt.plot(T410_P,x0,'k')
plt.plot(T420_P,x0,'y')
#plt.plot(T430_P,x0,'g')
plt.plot(T490_P,x0,'k')
#plt.plot(T520_P,x0,'b')
plt.plot(T660_P,x0,color='orange')
plt.plot(T1100_P,x0,color='gray')
#plt.plot(T425_P,x0,'y')
#plt.plot(T307_P,x0)
#plt.plot(T310_P,x0,'y')
#plt.plot(T313_P,x0)
#plt.plot(T320_P,x0,'k')
if show_arrows:
	for i in range(10,len(P0_323)):
		#ax.arrow(P0_323[i],R0_323[i],vector_323[0][i],vector_323[1][i],width=0.002,head_width=0.03,head_length=0.80,length_includes_head=True,fc='g',ec='g')
		#head_angle = math.atan2(vector_323[1][i],vector_323[0][i])
		#ax.arrow(P0_323[i],R0_323[i],vector_323[0][i],vector_323[1][i],head_length=abs(0.65*math.cos(head_angle)),head_width=abs(0.01*math.sin(head_angle)),length_includes_head=True,fc='g',ec='g')
		ax.arrow(P0_323[i],R0_323[i],vector_323[0][i],vector_323[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='g',ec='g')
		#ax.annotate("",xy=(1.5, 1.5),xytext=(0, 0),arrowprops=dict(arrowstyle="->"))
	for i in range(9,len(P0_373)):
		ax.arrow(P0_373[i],R0_373[i],vector_373[0][i],vector_373[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='r',ec='r')
	for i in range(0,len(P0_380)):
		ax.arrow(P0_380[i],R0_380[i],vector_380[0][i],vector_380[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='c',ec='c')
	for i in range(0,len(P0_400)):
		ax.arrow(P0_400[i],R0_400[i],vector_400[0][i],vector_400[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='m',ec='m')
	for i in range(0,len(P0_420)):
		ax.arrow(P0_420[i],R0_420[i],vector_420[0][i],vector_420[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='y',ec='y')
	for i in range(0,len(P0_490)):
		ax.arrow(P0_490[i],R0_490[i],vector_490[0][i],vector_490[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='k',ec='k')
	for i in range(0,len(P0_660)):
		ax.arrow(P0_660[i],R0_660[i],vector_660[0][i],vector_660[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='orange',ec='orange')
	for i in range(0,len(P0_1100)):
		ax.arrow(P0_1100[i],R0_1100[i],vector_1100[0][i],vector_1100[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='grey',ec='grey')
plt.xlabel('Pressure $P$ (MPa)',fontsize=axis_size)
plt.ylabel(r'Density $\rho$ (g/cm$^3$)',fontsize=axis_size)
legend = plt.legend(('Garg 323K','Garg 373K','Xiong 380K','Xiong 400K','Xiong 420K','Angus 490K','Angus 660K','Angus 1100K'),loc=2,fontsize=size,numpoints=1,ncol=2)
if not showbox:
	legend.get_frame().set_facecolor('none')
	legend.get_frame().set_linewidth(0.0)
plt.axis([0,70,0,1.2])
# minorLocator = AutoMinorLocator()
# ax.xaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
# minorLocator = AutoMinorLocator()
# ax.yaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
plt.tight_layout()
figPVSR.savefig('../'+output_folder+r'\pure_CO2_density_isotherm_'+param_set+img_extension,dpi=img_dpi,bbox_inches='tight')

#==================================================================================
# T versus R plots
figTVSR = plt.figure(num=None, figsize=(img_vert, img_hori), dpi=img_dpi, facecolor='w', edgecolor='k')
ax2 = plt.axes()
#plt.plot(T0_0pt5MPa,R0_0pt5MPa,'>k')
#plt.plot(T0_1MPa,R0_1MPa,'>k')
plt.plot(T0_6MPaG,R0_6MPaG,'or')
plt.plot(T0_12MPa,R0_12MPa,'sb')
plt.plot(T0_16MPa,R0_16MPa,'^g')
plt.plot(T0_30MPa,R0_30MPa,'vm')
plt.plot(T0_50MPa,R0_50MPa,'*c')
plt.plot(T0_65MPa,R0_65MPa,'<y')
#plt.plot(P0pt5MPa_T,y0,'k')
#plt.plot(P1MPa_T,y0,'k')
plt.plot(P6MPa_T,y0,'r')
plt.plot(P12MPa_T,y0,'b')
plt.plot(P16MPa_T,y0,'g')
plt.plot(P30MPa_T,y0,'m')
plt.plot(P50MPa_T,y0,'c')
plt.plot(P65MPa_T,y0,'y')
if show_arrows:
	for i in range(0,len(T0_30MPa),2):
		ax2.arrow(T0_30MPa[i],R0_30MPa[i],vector_30MPa[0][i],vector_30MPa[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='m',ec='m')
	for i in range(0,len(T0_50MPa),2):
		ax2.arrow(T0_50MPa[i],R0_50MPa[i],vector_50MPa[0][i],vector_50MPa[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='c',ec='c')
	for i in range(0,len(T0_65MPa),2):
		ax2.arrow(T0_65MPa[i],R0_65MPa[i],vector_65MPa[0][i],vector_65MPa[1][i],head_length=0,head_width=0,ls=arrow_ls,fc='y',ec='y')
plt.xlabel('Temperature $T$ (K)',fontsize=axis_size)
plt.ylabel(r'Density $\rho$ (g/cm$^3$)',fontsize=axis_size)
legend = plt.legend(('Vargaftik 6MPa','Vargaftik 12MPa','Vargaftik 16MPa','Vargaftik 30MPa','Angus 50MPa','Angus 65MPa'),loc=1,fontsize=size,numpoints=1)
if not showbox:
	legend.get_frame().set_facecolor('none')
	legend.get_frame().set_linewidth(0.0)
# minorLocator = AutoMinorLocator()
# ax2.xaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
# minorLocator = AutoMinorLocator()
# ax2.yaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
plt.axis([0,2000,0,1.20])
plt.tight_layout()
figTVSR.savefig('../'+output_folder+r'\pure_CO2_density_isobar_'+param_set+img_extension,dpi=img_dpi,bbox_inches='tight')

#==================================================================================
# R versus T coexistence plots.
figRVSTcoex = plt.figure(num=None, figsize=(img_vert, img_hori), dpi=img_dpi, facecolor='w', edgecolor='k')
ax3 = plt.axes()
plt.plot(T0_Coex_Newitt,R0_CoeG_Newitt,'^b')
plt.plot(T0_Coex_Angus,R0_CoeG_Angus,'ob')
plt.plot(T0_Coex_Duschek,R0_CoeG_Duschek,'sb')
plt.plot(T0_Coex_Newitt,R0_CoeL_Newitt,'^r')
plt.plot(T0_Coex_Angus,R0_CoeL_Angus,'or')
plt.plot(T0_Coex_Duschek,R0_CoeL_Duschek,'sr')
plt.plot(T0_Theory,rho_g,'-b')
plt.plot(T0_Theory,rho_l,'--r')
plt.xlabel('Temperature $T$ (K)',fontsize=axis_size)
plt.ylabel(r'Density $\rho$ (g/cm$^3$)',fontsize=axis_size)
legend = plt.legend(('Newitt Gas','Angus Gas','Duschek Gas','Newitt Liquid','Angus Liquid','Duschek Liquid'),loc=1,fontsize=size,numpoints=1,ncol=2)
if not showbox:
	legend.get_frame().set_facecolor('none')
	legend.get_frame().set_linewidth(0.0)
# minorLocator = AutoMinorLocator()
# ax3.xaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
# minorLocator = AutoMinorLocator()
# ax3.yaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
plt.axis([216.15,310,0,1.35])
plt.tight_layout()
figRVSTcoex.savefig('../'+output_folder+r'\pure_CO2_density_TR_coexistence_'+param_set+img_extension,dpi=img_dpi,bbox_inches='tight')

#==================================================================================
# R versus P coexistence plots.
figRVSP = plt.figure(num=None, figsize=(img_vert, img_hori), dpi=img_dpi, facecolor='w', edgecolor='k')
ax4 = plt.axes()
plt.plot(P0_Coex_Newitt,R0_CoeG_Newitt,'^b')
plt.plot(P0_Coex_Angus,R0_CoeG_Angus,'ob')
plt.plot(P0_Coex_Duschek,R0_CoeG_Duschek,'sb')
plt.plot(P0_Coex_Newitt,R0_CoeL_Newitt,'^r')
plt.plot(P0_Coex_Angus,R0_CoeL_Angus,'or')
plt.plot(P0_Coex_Duschek,R0_CoeL_Duschek,'sr')
plt.plot(P0_Theory,rho_g,'-b')
plt.plot(P0_Theory,rho_l,'--r')
plt.xlabel('Pressure $P$ (MPa)',fontsize=axis_size)
plt.ylabel(r'Density $\rho$ (g/cm$^3$)',fontsize=axis_size)
legend = plt.legend(('Newitt Gas','Angus Gas','Duschek Gas','Newitt Liquid','Angus Liquid','Duschek Liquid'),loc=1,fontsize=size,numpoints=1,ncol=2)
if not showbox:
	legend.get_frame().set_facecolor('none')
	legend.get_frame().set_linewidth(0.0)
# minorLocator = AutoMinorLocator()
# ax4.xaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
# minorLocator = AutoMinorLocator()
# ax4.yaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
plt.axis([0,8,0,1.35])
plt.tight_layout()
figRVSP.savefig('../'+output_folder+r'\pure_CO2_density_PR_coexistence_'+param_set+img_extension,dpi=img_dpi,bbox_inches='tight')

#==================================================================================
# Vapour pressure inversion curve plot.
T_Newitt_Inv = range(0,len(T0_Coex_Newitt))
P_Newitt_Inv = range(0,len(P0_Coex_Newitt))
for i in range(0,len(T0_Coex_Newitt)):
	T_Newitt_Inv[i] = 10**3/T0_Coex_Newitt[i]
	P_Newitt_Inv[i] = math.log(P0_Coex_Newitt[i])
T_Angus_Inv = range(0,len(T0_Coex_Angus))
P_Angus_Inv = range(0,len(P0_Coex_Angus))
for i in range(0,len(T0_Coex_Angus)):
	T_Angus_Inv[i] = 10**3/T0_Coex_Angus[i]
	P_Angus_Inv[i] = math.log(P0_Coex_Angus[i])
T_Duschek_Inv = range(0,len(T0_Coex_Duschek))
P_Duschek_Inv = range(0,len(P0_Coex_Duschek))
for i in range(0,len(T0_Coex_Duschek)):
	T_Duschek_Inv[i] = 10**3/T0_Coex_Duschek[i]
	P_Duschek_Inv[i] = math.log(P0_Coex_Duschek[i])
T_Theory_Inv = range(0,len(T0_Theory))
P_Theory_Inv = range(0,len(P0_Theory))
for i in range(0,len(T0_Theory)):
	T_Theory_Inv[i] = 10**3/T0_Theory[i]
	P_Theory_Inv[i] = math.log(P0_Theory[i])
figVAPINV = plt.figure(num=None, figsize=(img_vert, img_hori), dpi=img_dpi, facecolor='w', edgecolor='k')
ax5 = plt.axes()
plt.plot(T_Newitt_Inv,P_Newitt_Inv,'^g')
plt.plot(T_Angus_Inv,P_Angus_Inv,'om')
plt.plot(T_Duschek_Inv,P_Duschek_Inv,'sc')
plt.plot(T_Theory_Inv,P_Theory_Inv,'b')
plt.xlabel(r'$1/T$ (10$^{-3}$ K$^{-1}$)',fontsize=axis_size)
plt.ylabel(r'log$(P/P_0)$',fontsize=axis_size)
legend = plt.legend(('Newitt Gas','Angus Gas','Duschek Gas'),loc=1,fontsize=size,numpoints=1)
if not showbox:
	legend.get_frame().set_facecolor('none')
	legend.get_frame().set_linewidth(0.0)
# minorLocator = AutoMinorLocator()
# ax5.xaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
# minorLocator = AutoMinorLocator()
# ax5.yaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
plt.tight_layout()
figVAPINV.savefig('../'+output_folder+r'\pure_CO2_vapour_pressure_inversion_'+param_set+img_extension,dpi=img_dpi,bbox_inches='tight')

#==================================================================================
# Thermal expansivity plots.
figTHREXP = plt.figure(num=None, figsize=(img_vert, img_hori), dpi=img_dpi, facecolor='w', edgecolor='k')
ax6 = plt.axes()
plt.plot(T0_16MPa,A0_16MPa,'or')
plt.plot(T0_30MPa,A0_30MPa,'sg')
plt.plot(T0_50MPa,A0_50MPa,'db')
plt.plot(T0_65MPa,A0_65MPa,'*y')
plt.plot(P16MPa_T,P16MPa_A,'r')
plt.plot(P30MPa_T,P30MPa_A,'g')
plt.plot(P50MPa_T,P50MPa_A,'b')
plt.plot(P65MPa_T,P65MPa_A,'y')
plt.xlabel('Temperature $T$ (K)',fontsize=axis_size)
plt.ylabel(r'Thermal expansivity $\alpha_V$ (K$^{-1}$)',fontsize=axis_size)
legend = plt.legend(('Vargaftik 16MPa','Vargaftik 30MPa','Angus 50MPa','Angus 65MPa'),loc=1,fontsize=size,numpoints=1)
if not showbox:
	legend.get_frame().set_facecolor('none')
	legend.get_frame().set_linewidth(0.0)
plt.axis([210,1200,0,0.018])
# minorLocator = AutoMinorLocator()
# ax6.xaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
# minorLocator = AutoMinorLocator()
# ax6.yaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
plt.tight_layout()
figTHREXP.savefig('../'+output_folder+r'\pure_CO2_thermal_expansivity_'+param_set+img_extension,dpi=img_dpi,bbox_inches='tight')

#==================================================================================
#Isothermal compressibility plots.
figTHRCOMP = plt.figure(num=None, figsize=(img_vert, img_hori), dpi=img_dpi, facecolor='w', edgecolor='k')
ax7 = plt.axes()
plt.plot(P0_380,B0_380,'or')
plt.plot(P0_380A,B0_380A,'sr')
plt.plot(P0_430,B0_430,'dg')
plt.plot(P0_490,B0_490,'*b')
plt.plot(P0_660,B0_660,'<m')
plt.plot(P0_1100,B0_1100,'^c')
plt.plot(T380_P,T380_B,'r')
plt.plot(T430_P,T430_B,'g')
plt.plot(T490_P,T490_B,'b')
plt.plot(T660_P,T660_B,'m')
plt.plot(T1100_P,T1100_B,'c')
plt.xlabel('Pressure $P$ (MPa)',fontsize=axis_size)
plt.ylabel(r'Isothermal compressibility $\beta_T$ (MPa$^{-1}$)',fontsize=axis_size)
legend = plt.legend(('Xiong 380K','Angus 380K','Angus 430K','Angus 490K','Angus 660K','Angus 1100K'),loc=1,fontsize=size,numpoints=1)
if not showbox:
	legend.get_frame().set_facecolor('none')
	legend.get_frame().set_linewidth(0.0)
plt.axis([0,68,0,0.4])
# minorLocator = AutoMinorLocator()
# ax7.xaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
# minorLocator = AutoMinorLocator()
# ax7.yaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
plt.tight_layout()
figTHRCOMP.savefig('../'+output_folder+r'\pure_CO2_isothermal_compressibility_'+param_set+img_extension,dpi=img_dpi,bbox_inches='tight')

#==================================================================================
# H_vap versus T coexistence plots.
figHVAP = plt.figure(num=None, figsize=(img_vert, img_hori), dpi=img_dpi, facecolor='w', edgecolor='k')
ax10 = plt.axes()
plt.plot(T0_Coex_Newitt,H0_Coex_Newitt,'^r')
plt.plot(T0_Coex_Angus,H0_Coex_Angus,'og')
plt.plot(T0_Theory,H_vap,'b')
plt.xlabel('Temperature $T$ (K)',fontsize=axis_size)
plt.ylabel(r'Enthalpy of vapourization $\Delta H_{vap}$ (J/g)',fontsize=axis_size)
legend = plt.legend(('Newitt','Angus'),loc=1,fontsize=size,numpoints=1)
if not showbox:
	legend.get_frame().set_facecolor('none')
	legend.get_frame().set_linewidth(0.0)
#plt.title(title,fontsize=title_size)
#plt.axis([216.15,310,0,1.2])
# minorLocator = AutoMinorLocator()
# ax10.xaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
# minorLocator = AutoMinorLocator()
# ax10.yaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
plt.tight_layout()
figHVAP.savefig('../'+output_folder+r'\pure_CO2_enthalpy_vapourization_'+param_set+img_extension,dpi=img_dpi,bbox_inches='tight')

#==================================================================================
# Second virial coefficient along zero pressure isobar plot.
figVIRIAL = plt.figure(num=None, figsize=(img_vert, img_hori), dpi=img_dpi, facecolor='w', edgecolor='k')
ax11 = plt.axes()
plt.plot(T0_0MPa,BC_0MPa,'ob')
plt.plot(T_SVC,P0MPa_BC,'g')
plt.xlabel('Temperature $T$ (K)',fontsize=axis_size)
plt.ylabel(r'Second virial coefficient $B$ (cm$^3$/g)',fontsize=axis_size)
#plt.legend(('Angus 0MPa','Theory'),loc=4,fontsize=size,numpoints=1)
plt.axis([0,1200,-6,1])
# minorLocator = AutoMinorLocator()
# ax11.xaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
# minorLocator = AutoMinorLocator()
# ax11.yaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
plt.tight_layout()
figVIRIAL.savefig('../'+output_folder+r'\pure_CO2_second_virial_coefficient_'+param_set+img_extension,dpi=img_dpi,bbox_inches='tight')

#==================================================================================
#Relative deviation at single temperature plot.
figTREL2 = plt.figure(num=None, figsize=(img_vert_dev, img_hori_dev), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()
exec "plt.plot(P0_%s,Kilpatrick_deviation,'-^r')" % (temp)
exec "plt.plot(P0_%s,Kiszka_deviation,'-vg')" % (temp)
exec "plt.plot(P0_%s,Pope_deviation,'-xm')" % (temp)
exec "plt.plot(P0_%s,Wang_deviation,'-d',color='Sienna')" % (temp)
exec "plt.plot(P0_%s,Hariharan_deviation,'-sy')" % (temp)
exec "plt.plot(P0_%s,Garg_deviation,'-vk')" % (temp)
exec "plt.plot(P0_%s,Xiong_deviation,'-^',color='GreenYellow')" % (temp)
exec "plt.plot(P0_%s,Doghieri_deviation,'-x',color='gray')" % (temp)
exec "plt.plot(P0_%s,Nalawade_deviation,'-p',color='Orange')" % (temp)
exec "plt.plot(P0_%s,Funami_deviation,'-<',color='DarkOliveGreen')" % (temp)
exec "plt.plot(P0_%s,Cao_deviation,'->',color='FireBrick')" % (temp)
exec "plt.plot(P0_%s,Arce_deviation,'-*c')" % (temp)
exec "plt.plot(P0_%s,Self_deviation,'-ob')" % (temp)

plt.xlabel(r'Pressure $P$ (MPa)',fontsize=axis_size)
plt.ylabel(r'Relative density deviation (%)',fontsize=axis_size)
legend = plt.legend(('Kilpatrick','Kiszka','Pope','Wang','Hariharan','Garg','Xiong','Doghieri','Nalawade','Funami','Cao','Arce','Present work'),loc=2,fontsize=size,numpoints=1,ncol=4)
if not showbox:
	legend.get_frame().set_facecolor('none')
	legend.get_frame().set_linewidth(0.0)
plt.axis([0,70,0,30])
# minorLocator = AutoMinorLocator()
# ax.xaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
# minorLocator = AutoMinorLocator()
# ax.yaxis.set_minor_locator(minorLocator)
# plt.tick_params(which='both', width=1)
# plt.tick_params(which='major', length=7)
# plt.tick_params(which='minor', length=4)
plt.tight_layout()
figTREL2.savefig('../'+output_folder+r'\pure_CO2_temperature_relative_deviation_single'+img_extension,dpi=img_dpi,bbox_inches='tight')
