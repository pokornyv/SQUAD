# SQUAD - superconducting quantum dot
# reading and parsing parameter file ssn.in
# Vladislav Pokorny, 2016; pokornyv@fzu.cz

import scipy as sp
from ConfigParser import SafeConfigParser

config = SafeConfigParser()
config.optionxform = str		# case-sensitive names
config.read('ssn.in')

P = {}

# default values ##########################################

P['M']            = 21
P['dE']           = 1e-4
P['rootf']        = 'brentq'
P['ConvN']        = 1e-4
P['ConvX']        = 1e-5
P['ConvHF']       = 1e-6
P['MuMin']        = -2.0
P['MuMax']        = 2.0
P['HF_max_iter']  = 10000
P['offset_x']     = 1e-12

P['WriteIO']      = 1
P['Write_HFGF']   = 0
P['Write_Bubble'] = 0
P['Write_2ndSE']  = 0
P['Write_2ndGF']  = 1
P['Write_AC']     = 0
P['EmaxFiles']    = 10.0
P['EstepFiles']   = 10

# read in file ############################################

params_int_L = ['M','HF_max_iter','EstepFiles']
params_float_L = ['dE','ConvN','ConvX','ConvHF','MuMin','MuMax','offset_x','EmaxFiles']
params_bool_L = ['WriteIO','Write_HFGF','Write_Bubble','Write_2ndSE','Write_2ndGF','Write_AC']

for sec in config.sections():
	for opt in config.options(sec):
		if opt in params_int_L:
			P[opt] = int(config.get(sec,opt))
		elif opt in params_float_L:
			P[opt] = float(config.get(sec,opt))
		elif opt in params_bool_L:
			P[opt] = bool(int(config.get(sec,opt)))
		else:
			P[opt] = str(config.get(sec,opt))

## params.py end ##

