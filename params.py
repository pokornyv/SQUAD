import scipy as sp
from ConfigParser import SafeConfigParser

config = SafeConfigParser()
config.read('ssn.in')

M  = int(config.get('params','M'))
dE = float(config.get('params','dE'))

rootf = str(config.get('params','rootf'))

ConvN       = float(config.get('params','ConvN'))
ConvX       = float(config.get('params','ConvX'))
ConvHF      = float(config.get('params','ConvHF'))
MuMin       = float(config.get('params','MuMin'))
MuMax       = float(config.get('params','MuMax'))
HF_max_iter = int(config.get('params','HF_max_iter'))

chat       = bool(int(config.get('params','WriteIO')))
WriteFiles = bool(int(config.get('params','WriteFiles')))
EmaxFiles  = float(config.get('params','EmaxFiles'))
EstepFiles = int(config.get('params','EstepFiles'))

offset_x  = float(config.get('params','offset_x'))

## params.py end ##

