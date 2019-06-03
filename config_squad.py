################################################################
# SQUAD - superconducting quantum dot                          #
# Copyright (C) 2012-2019  Vladislav Pokorny; pokornyv@fzu.cz  #
# homepage: github.com/pokornyv/SQUAD                          #
# config_squad.py - config and global variables                #
# method described in                                          #
#    Sci. Rep. 5, 8821 (2015).                                 #
#    Phys. Rev. B 93, 024523 (2016).                           #
################################################################

from __future__ import print_function
import scipy as sp
from os import listdir
from sys import argv,exit,version_info
from time import ctime,time
from ConfigParser import SafeConfigParser

###########################################################
## reading parameters from command line ###################

U      = float(argv[1])
Delta  = float(argv[2])
GammaR = float(argv[3])
GammaL = float(argv[4])*GammaR
eps    = float(argv[5])
P      = float(argv[6])
GammaN = 0.0 	 ## for compatibility with functions from SSN branch (not in GitHub master branch)

## little post-processing
ed       = eps-U/2.0                               ## energy level shifted to symmetry point
Phi      = P*sp.pi                                 ## phase difference
GammaLR  = GammaL/GammaR if GammaR != 0.0 else 1.0 ## coupling strength ratio
GammaTot = GammaL + GammaR                         ## total coupling strength

###########################################################
## reading config file ####################################

cfile = 'squad.in'

if cfile not in listdir('.'): 
	print('- Parameter file '+cfile+' missing. Exit.')
	exit(1)

config = SafeConfigParser()
config.optionxform = str    ## case-sensitive names
config.read(cfile)

## default values #########################################
## consult infile.md for details

M            = 20
dE           = 1e-4
rootf        = 'brentq'
ConvN        = 1e-4
ConvX        = 1e-5
ConvHF       = 1e-6
MuMin        = -2.0
MuMax        = 2.0
ABSinit_val  = 0.99
HF_max_iter  = 10000
offset_x     = 1e-12

WriteIO      = 1
Write_HFGF   = 0
Write_Bubble = 0
Write_2ndSE  = 0
Write_2ndGF  = 1
Write_AC     = 0           ## Andreev conductance, for compatibility with SSN codes
EmaxFiles    = 10.0
EstepFiles   = 10

## read the .in file ######################################

## [params] section
if config.has_option('params','M'):
	M            = int(config.get('params','M'))
if config.has_option('params','dE'):
	dE           = float(config.get('params','dE'))
if config.has_option('params','rootf'):
	rootf        = str(config.get('params','rootf'))
if config.has_option('params','ConvN'):
	ConvN        = float(config.get('params','ConvN'))
if config.has_option('params','ConvX'):
	ConvX        = float(config.get('params','ConvX'))
if config.has_option('params','ConvHF'):
	ConvHF       = float(config.get('params','ConvHF'))
if config.has_option('params','MuMin'):
	MuMin        = float(config.get('params','MuMin'))
if config.has_option('params','MuMax'):
	MuMax        = float(config.get('params','MuMax'))
if config.has_option('params','ABSinit_val'):
	ABSinit_val  = float(config.get('params','ABSinit_val'))
if config.has_option('params','HF_max_iter'):
	HF_max_iter  = int(config.get('params','HF_max_iter'))
if config.has_option('params','offset_x'):
	offset_x     = float(config.get('params','offset_x'))

## [IO] section
if config.has_option('IO','WriteIO'):
	chat         = bool(int(config.get('IO','WriteIO')))
if config.has_option('IO','Write_HFGF'):
	Write_HFGF   = bool(int(config.get('IO','Write_HFGF')))
if config.has_option('IO','Write_Bubble'):
	Write_Bubble = bool(int(config.get('IO','Write_Bubble')))
if config.has_option('IO','Write_2ndSE'):
	Write_2ndSE  = bool(int(config.get('IO','Write_2ndSE')))
if config.has_option('IO','Write_2ndGF'):
	Write_2ndGF  = bool(int(config.get('IO','Write_2ndGF')))
if config.has_option('IO','Write_AC'):
	Write_AC     = bool(int(config.get('IO','Write_AC')))		## compatibility with SSN codes
if config.has_option('IO','EmaxFiles'):
	EmaxFiles    = float(config.get('IO','EmaxFiles'))
if config.has_option('IO','EstepFiles'):
	EstepFiles   = int(config.get('IO','EstepFiles'))

###########################################################
## energy axis ############################################

## In case you run into RuntimeWarning: invalid value encountered in power:
## for Kramers-Kronig we need range(N)**3 array, for large N it can 
## hit the limit of 2**63 = 9223372036854775808 of signed int
## large values of N also introduce instability to calcualtion of ABS
N      = 2**M-1	## number of points for bubble/self-energy FFT calculation
dE_dec = int(-sp.log10(dE))
En_A   = sp.around(sp.linspace(-(N-1)/2*dE,(N-1)/2*dE,N),dE_dec+2)
Nhalf  = int((len(En_A)-1)/2)	## zero on the energy axis

## cannot print what is not calculated
if EmaxFiles > sp.fabs(En_A[0]): EmaxFiles = sp.fabs(En_A[0])  

if any([GammaL <= 0.0,GammaR <= 0.0, U < 0.0, Delta <= 0.0]):
	print('# check_params: Error: All of GammaL, GammaR, U, Delta must be positive.')
	exit(1)

if Delta > En_A[-1]:
	print('# Error: Delta must be smaller than the bandwidth.')
	print('# Delta = {0: .5f}, Emax = {1: .5f}'.format(Delta,En_A[-1]))
	exit(1)

## locate band edges in En_A
EdgePos1 = sp.nonzero(sp.around(En_A,dE_dec) == sp.around(-Delta,dE_dec))[0][0]
EdgePos2 = sp.nonzero(sp.around(En_A,dE_dec) == sp.around( Delta,dE_dec))[0][0]

## Fermi-Dirac distribution for T=0
FD_A = 1.0*sp.concatenate([sp.ones(int((N-1)/2)),[0.5],sp.zeros(int((N-1)/2))])

## config_squad.py end ##

