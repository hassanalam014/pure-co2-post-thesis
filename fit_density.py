# Date: April 2017
#
# Description: The purpose of this file is to estimate the pure fluid parameters of Carbon Dioxide (CO2) by regressing to experimental PVT data.
#

import os,sys,math,numpy as npy
from loadExperimentalData import *
from lmfit import minimize, Parameters, report_fit
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
from calculatePureResidual import calculatePureResidual,excludeCP

N = len(P0)
deltaP = max(P0)-min(P0)
deltaT = max(T0)-min(T0)
print('Performing fit with {} datapoints over a temperature range of {}-{}K and a pressure range of {}-{}MPa.'.format(N,min(T0),max(T0),round(min(P0),2),max(P0)))

PE = 3.0
TE = 15.0

exclCP = True

P0,T0,R0,M0,I0 = excludeCP(exclCP,P0,T0,R0,M0,I0,Pc0,Tc0,Rc0,PE,TE)

#Initializing the parameters.
params = Parameters()
#The following code sets up the model's parameters. It includes both fitting parameters and parameters that will remain fixed
#for the fitting. The values given are the inital guesses of fitting parameters and values of fixed parameters.
#				(Name,		Value,		Vary?,	Min,	Max,	Expr)
params.add_many(('Pstar',	419.9,		True,	0,		None,	None),
				('Tstar',	341.8,		True,	0,		None,	None),
				('Rstar',	1.397,		True,	0,		None,	None),
				('Pc0',		Pc0,		False,	0,		None,	None),
				('Tc0',		Tc0,		False,	0,		None,	None),
				('Rc0',		Rc0,		False,	0,		None,	None))

#Running the Levenberg-Marquart algorithm on the residuals in order to do least squares fitting. This will return the fitted value of the RESIDUALS.
#These need to be added to the experimental datapints to find the fitted pressures.
fit = minimize(calculatePureResidual,params,args=(P0,T0,R0,M0,I0,PE,TE,'P'))

#Reporting the values of the parameters. NEED TO FIGURE OUT HOW TO PRINT THIS TO FILE.
report_fit(fit.params)
