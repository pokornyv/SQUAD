# script for running codes for a selected set of input parameters

from os import system
from sys import argv,exit
import scipy as sp

U       = 1.0*eval(argv[1])
Delta   = 1.0*eval(argv[2])
GammaR  = 1.0*eval(argv[3])
GammaLR = 1.0*eval(argv[4])
GammaN  = 1.0*eval(argv[5])
en      = 1.0*eval(argv[6])
P       = 1.0*eval(argv[7])

xmin = 1.0*eval(argv[8])
xmax = 1.0*eval(argv[9])
dx   = 1.0*eval(argv[10])

# choose the parameter: 1- U, 2 - E, 3 - P, 4 - GammaR - interaction, local energy, phase difference, right coupling
# the previous set of this parameter acts as dummy 
par = eval(argv[11])

x = xmin
while x < xmax + dx:
	if   par == 1: U      = x
	elif par == 2: en     = x
	elif par == 3: P      = x
	elif par == 4: GammaR = x
	else: GammaN = x
	comm="python ssn_second.py "+str(U)+" "+str(Delta)+" "+str(GammaR)+" "+str(GammaLR)+" "+str(GammaN)+" "+str(en)+" "+str(P)
	k = system(comm)
	if k != 0: break
	x = x + dx

