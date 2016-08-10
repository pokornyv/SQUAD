# script for running codes for a selected set of input parameters

from os import system
from sys import argv,exit
import scipy as sp

U       = 1.0*eval(argv[1])
Delta   = 1.0*eval(argv[2])
GammaR  = 1.0*eval(argv[3])
GammaLR = 1.0*eval(argv[4])
en      = 1.0*eval(argv[5])
P       = 1.0*eval(argv[6])

xmin = 1.0*eval(argv[7])
xmax = 1.0*eval(argv[8])
dx   = 1.0*eval(argv[9])

# choose the parameter: 1- U, 2 - E, 3 - P, 4 - GammaR - interaction, local energy, phase difference, right coupling
# the previous set of this parameter acts as dummy 
par = eval(argv[10])

x = xmin
while x < xmax + dx:
	if   par == 1: U      = x
	elif par == 2: en     = x
	elif par == 3: P      = x
	elif par == 4: GammaR = x
	else: Gamma = x
	comm="python second_fullsc.py "+str(U)+" "+str(Delta)+" "+str(GammaR)+" "+str(GammaLR)+" "+str(en)+" "+str(P)+" 0"
	k = system(comm)
	if k != 0: break
	x = x + dx

