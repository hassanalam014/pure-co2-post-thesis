# Date: July 2017
#
# Description: The purpose of this file is caulculate the free energy and equation of state of Carbon Dioxide
# 			   based on the regression performed in the paper:
#			   Span et al., 1996. Journal of Physical and Chemical Reference Data (25).
#			   #Valid from the triple point 216.592K to 1100K and up to 800MPa.
#			   Not valid near the critical point (7.377MPa, 304.128K).
#

import os,sys,math,numpy as npy
from sympy import *

#
# isListOrNpyArray( variable )
#
# Summary		: The function determines the shape of an
#				  input variable to determine if it is vector or scalar.
#
# Parameters	: variable: an integer, float, npy.float, array, npy.array
#
# Returns		: iLONA: Boolean (True/False)
#

def isListOrNpyArray(variable):

	if isinstance(variable,list) or isinstance(variable,npy.ndarray):
		iLONA = True
	else:
		iLONA = False
		
	return iLONA

#
# HelmholtzCO2( T, R )
#
# Summary		: The function calculates the Helmholtz free energy per
#				  molecule F/N of Carbon Dioxide based on the free energy function
#				  given by Span et al, 1996. Journal of Physical and Chemical Reference Data (25).
#
# Parameters	: T: a float representing temperature.
#				: R: a float representing density.
#
# Returns		: phi_value: a float representing free energy per molecule.
#

def HelmholtzDME(T,R):
	#Molar mass of Carbon Dioxide.
	M = 44.01
	#Critical temperature in K.
	T_c = 304.128
	#Critical density.
	R_c = 0.4676

	print R

	#Reduced temperature and density parameters.
	tau0 = T_c/T
	delta0 = R/R_c

	#Setting up variable for Sympy.
	tau = Symbol('tau')
	delta = Symbol('delta')

	#IDEAL GAS FREE ENERGY COMPONENT.
	#Ideal gas model parameters (corresponding to phi^0).
	a0 = [8.37304456,-3.70454304,2.50000000,1.99427042,0.62105248,0.41195293,1.04028922,0.08327678]
	theta0 = [0,0,0,3.15163,6.11190,6.77708,11.32384,27.08792]

	#Ideal gas reduced free energy.
	phi0 = log(delta)+a0[0]+a0[1]*tau+a0[2]*log(tau)
	for i in range(3,8):
		phi0 += a0[i]*log(1.0-exp(-tau*theta0[i]))

	#RESIDUAL FREE ENERGY COMPONENT.

	#Residual parameters of the nonanalytic terms.
	ai = [3.00,3.50,4.00]
	bj = [0.875,0.925,0]
	Bk = [0.30,1.00,0]
	Cl = [10.00,12.50,15.00]
	Dm = [225.0,250.0,275.0]
	A = 0.700
	beta = 0.300

	#Residual model parameters (corresponding to phi^r).
	n = [0.38856823203161,2.9385475942740,-5.5867188534934,-0.76753199592477,0.31729005580416,0.54803315897767,0.12279411220335,2.1658961543220,1.5841735109724,-0.23132705405503,
			0.058116916431436,-0.55369137205382,0.48946615909422,-0.024275739843501,0.062494790501678,-0.12175860225246,-0.37055685270086,-0.016775879700426,-0.11960736637987,-0.045619362508778,
			0.035612789270346,-0.0074427727132052,-0.0017395704902432,-0.021810121289527,]
	d = []
	t = []
	c = []
	b = []
	g = []
	e = []

	#Variables used to 'toggle' each of the three terms in the sum.
	term1 = [1,1,1,1,1,0,0,0,0,0,0,0,0,0,0]
	term2 = [0,0,0,0,0,1,1,1,1,1,1,0,0,0,0]
	term3 = [0,0,0,0,0,0,0,0,0,0,0,1,1,1,1]

	#Ideal gas reduced free energy.
	phir = 0
	for i in range(0,15):
		phir += term1[i]*n[i]*delta**d[i]*tau**t[i]
		phir += term2[i]*n[i]*delta**d[i]*tau**t[i]*exp(-delta**l[i])
		phir += term3[i]*n[i]*delta**d[i]*tau**t[i]*exp(-h[i]*(delta-e[i])**2-b[i]*(tau-g[i])**2)

	phi = phi0 + phir

	phi_value = phi.subs([(delta,delta0),(tau,tau0)])

	return phi_value

#
# pressureDME( T, R )
#
# Summary		: The function calculates the pressure based on the equation
#				  of state derived from Wu et al, 2011. Journal of Physical and Chemical Reference Data (40).
#
# Parameters	: T: a float representing temperature.
#				: R: a float representing density.
#
# Returns		: P: a float representing pressure.
#

def pressureDME(T,R):
	#Molar mass of Carbon Dioxide.
	M = 44.01
	#Critical temperature in K.
	T_c = 400.378
	#Critical density.
	R_c = 5.94

	#Reduced temperature and density parameters.
	tau0 = T_c/T
	delta0 = R/R_c

	#Setting up variable for Sympy.
	tau = Symbol('tau')
	delta = Symbol('delta')

	#IDEAL GAS FREE ENERGY COMPONENT.
	#Ideal gas model parameters (corresponding to phi^0).
	a = [-1.980976,3.171218]
	c0 = 4.039
	v = [2.641, 2.123, 8.992, 6.191]
	u = [361.0, 974.0, 1916.0, 4150.0]

	#Ideal gas reduced free energy.
	phi0 = a[0]+a[1]*tau+log(delta)+(c0-1.0)*log(tau)
	for i in range(0,4):
		phi0 -= v[i]*log(1.0-exp(-u[i]*tau/T_c))

	#RESIDUAL FREE ENERGY COMPONENT.
	#Residual model parameters (corresponding to phi^r).
	n = [0.029814139, 1.43517, -2.64964, -0.29515532, 0.17035607, -0.94642918, -0.099250514, 1.1264071, -0.76936548, -0.020717696, 0.24527037, 1.1863438, -0.49398368, -0.16388716, -0.027583584]
	d = [4.0, 1.0, 1.0, 2.0, 3.0, 1.0, 3.0, 2.0, 2.0, 7.0, 1.0, 1.0, 1.0, 3.0, 3.0]
	t = [1.0, 0.4366, 1.011, 1.137, 0.45, 2.83, 1.5, 1.235, 2.675, 0.7272, 1.816, 1.783, 3.779, 3.282, 1.059]
	l = [0, 0, 0, 0, 0, 2.0, 2.0, 1.0, 2.0, 1.0, 1.0, 0, 0, 0, 0]
	h = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.965336, 1.50858, 0.963855, 9.72643]
	b = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.28719, 0.806235, 0.777942, 197.681]
	g = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.277720, 0.430750, 0.429607, 1.138490]
	e = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.672698, 0.924246, 0.750815, 0.800022]

	#Variables used to 'toggle' each of the three terms in the sum.
	term1 = [1,1,1,1,1,0,0,0,0,0,0,0,0,0,0]
	term2 = [0,0,0,0,0,1,1,1,1,1,1,0,0,0,0]
	term3 = [0,0,0,0,0,0,0,0,0,0,0,1,1,1,1]

	#Ideal gas reduced free energy.
	phir = 0
	for i in range(0,15):
		phir += term1[i]*n[i]*delta**d[i]*tau**t[i]
		phir += term2[i]*n[i]*delta**d[i]*tau**t[i]*exp(-delta**l[i])
		phir += term3[i]*n[i]*delta**d[i]*tau**t[i]*exp(-h[i]*(delta-e[i])**2-b[i]*(tau-g[i])**2)

	phi = phi0 + phir

	pressure = 8.314472*T*R_c*delta**2*diff(phi,delta)*delta

	P = pressure.subs([(delta,delta0),(tau,tau0)])

	return P

#
# calculatePressureDME( T0, R0 )
#
# Summary		: A wrapper function for pressureDME that customizes the shape of the output
#				  pressure based on the shape of the input arrays.
#
# Parameters	: T0: a float or array representing temperature.
#				: R0: a float or array representing density.
#
# Returns		: P: a float or array representing pressure.
#

def calculatePressureDME(T0,R0):

    if not isListOrNpyArray(T0) and not isListOrNpyArray(R0):
        
        P = pressureDME(T0,R0)

        result = P
    elif isListOrNpyArray(R0):
        P = range(0,len(R0))

        for i in range(0,len(R0)):
            P[i] = pressureDME(T0,R0[i])

        result = P
    elif isListOrNpyArray(T0):  
        P = range(0,len(T0))

        for i in range(0,len(T0)):
            P[i] = pressureDME(T0[i],R0)

        result = P

    return result
