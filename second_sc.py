# SQUAD - superconducting quantum dot
# self-consistent second order perturbation theory solver
# for general case of asymmetric couplings GammaL != GammaR
# uses scipy, optimized on Python 2.7.5
# Vladislav Pokorny, 2015; pokornyv@fzu.cz

#header for output files:
#U		GammaR	GammaL	eps		Phi/pi	wABS		n			mu			ResGn1		ResGn2		ResA1

import scipy as sp
from scipy.optimize import brentq,fixed_point
from sys import argv,exit,version_info
from time import ctime
from squadlib1 import FillEnergies2,AndreevEnergy,GFresidues,SolveHF,FillGreenHF
from squadlib2 import *

U      = 1.0*eval(argv[1])
Delta  = 1.0*eval(argv[2])
GammaR = 1.0*eval(argv[3])
GammaL = 1.0*eval(argv[4])*GammaR
eps    = 1.0*eval(argv[5])
P      = 1.0*eval(argv[6])
try: chat  = eval(argv[7])
except IndexError: chat = 1	# 1 = write progress lines to standard output

#print '{0: .14f} {1: .14f}'.format(GammaR,GammaL)

WriteOut = 0				# 1 = write datafiles
ed  = eps-U/2.0				# localized energy level shifted to symmetry point
Phi = P*sp.pi
Conv = 2e-4					# convergence criterium for n and mu
ConvX = 1e-5				# convergence criterium for brentq/fixed_point

FitMin = 20.0
FitMax = 30.0

#N  = 2**19-1				# number of points for bubble/self-energy fft calculation
#dE = 1e-4

# RuntimeWarning: invalid value encountered in power:
# for KK we need range(N)**3, for large arrays it can 
# hit the limit of 9223372036854775808 == 2**63 of signed int
# large values of N also introduce instability to calcualtion of ABS

N  = 2**19-1				# number of points for bubble/self-energy fft calculation
dE = 1e-4
dE_dec = int(-sp.log10(dE))
En_F   = FillEnergies2(dE,N)
EminExp = 10			# maximum of the frequency range for outputs
SEtype  = 'sc2nd'		# identifier for output files

#rootf = 'fixed_point'	# root finding algorithm for density iterations # fixed point is not tested!!!
rootf = 'brentq'		# starts from wrong first guess

# calculating the Hartree-Fock parameters #################
if chat: print '#######################################################################################'
ver = str(version_info[0])+'.'+str(version_info[1])+'.'+str(version_info[2])
if chat: print '# generated by '+str(argv[0])+', python version: '+str(ver)+\
', SciPy version: '+str(sp.version.version)+', '+str(ctime())
if chat: print '# U ={0: .3f}, Delta ={1: .3f}, GammaR ={2: .3f}, GammaL ={3: .3f}, eps ={4: .3f}, Phi/pi ={5: .3f}'\
.format(U,Delta,GammaR,GammaL,eps,P)

if chat: print '# calculating HF solution...'
try:
	[n,mu,wzero,ErrMsgHF] = SolveHF(U,Delta,GammaR,GammaL,eps,P)
except RuntimeError:
	print '#  Warning: failed to calculate HF solution.'
	exit(0)

params_F = [U,Delta,GammaR,GammaL,P,eps]
hfe = ed+U*n				# Hartree-Fock energy
wzero = AndreevEnergy(U,GammaR,GammaL,Delta,Phi,hfe,mu)		# HF ABS frequencies
#print '{0: .3f}\t{1: .3f}\t{2: .3f}\t{3: .3f}\t{4: .3f}\t{5: .5f}\t{6: .5f}\t{7: .5f}'\
#.format(U,GammaR,GammaL,eps,P,float(wzero),float(n),float(mu)) # HF solution

if chat: print '# initial HFA values: n ={0: .5f}, mu ={1: .5f}, wzero ={2: .5f}'\
.format(float(n),float(mu),float(wzero))

[GFn_F,GFa_F,EdgePos1,EdgePos2,ABSposGF1,ABSposGF2] = FillGreenHF(U,Delta,GammaR,GammaL,hfe,Phi,mu,wzero,En_F)
#if WriteOut == 1: WriteFile(En_F,GFn_F,GFa_F,params_F,wzero,'HF_green',EminExp,10)
[ResGnp1,ResGnh1,ResGa1] = GFresidues(U,Delta,GammaR,GammaL,hfe,Phi,mu,-wzero)	# HF residues
[ResGnp2,ResGnh2,ResGa2] = GFresidues(U,Delta,GammaR,GammaL,hfe,Phi,mu, wzero)
ResOld_F = sp.array([ResGnp1,ResGnp2,ResGa1])
IDin = IntDOS(GFn_F,En_F)

if chat: print '# - energy interval minimum: '+str(En_F[0])
if chat: print '# - input from HF: w(ABS) ={0: .5f}, intDOS ={1: .5f}'.format(float(wzero),float(IDin))
if chat: print '# - residues: Gn: [{0: .5f}, {1: .5f}], Ga: {2: .5f}'.format(float(ResGnp1),float(ResGnp2),float(ResGa1))
if WriteOut == 1: WriteFile(En_F,GFn_F,GFa_F,params_F,wzero,'HF_green',EminExp,10)

# bubbles and vertex ######################################
[Chin_F,Chia_F,ABSposChi1,ABSposChi2] = TwoParticleBubbles(GFn_F,GFa_F,En_F,wzero)
Chin_F = FitTail(En_F,Chin_F,FitMin,FitMax,'even')
if WriteOut == 1: WriteFile(En_F,Chin_F,Chia_F,params_F,En_F[ABSposChi1],'HF_bubbles',EminExp,10)
if SEtype == 'sc2nd': ChiGamma_F = U**2*(Chin_F+Chia_F)  # second-order kernel of SDE (without HF)

# dynamical self-energy ###################################
[Sigman_F,Sigmaa_F] = SelfEnergy(GFn_F,GFa_F,ChiGamma_F,En_F)
Sigman_F = FitTail(En_F,Sigman_F,FitMin,FitMax,'odd')
Sigmaa_F = FitTail(En_F,Sigmaa_F,FitMin,FitMax,'even')
if WriteOut == 1: WriteFile(En_F,Sigman_F,Sigmaa_F,params_F,0.0,'2nd_SE',EminExp,10)

# initial guess for static part of self-energy ############
n  = ElectronDensity(params_F,n,mu,Sigman_F,Sigmaa_F,En_F)
mu = CooperPairDensity(params_F,n,mu,Sigman_F,Sigmaa_F,En_F)

# static self-energy ######################################
if chat: print '# - iterating densities'
n_old = 1e5
mu_old = 1e5
while any([sp.fabs(n-n_old)>Conv,sp.fabs(mu-mu_old)>Conv]):
	n_old = n
	mu_old = mu
	if rootf == 'brentq':
		if eps == 0.0: n = 0.5
		else: 
			eqn1 = lambda x: x - ElectronDensity(params_F,x,mu,Sigman_F,Sigmaa_F,En_F)
			n = brentq(eqn1,0.0,1.0,xtol=ConvX)
		eqn2 = lambda x: x - CooperPairDensity(params_F,n,x,Sigman_F,Sigmaa_F,En_F)
		mu = brentq(eqn2,-2.0,2.0,xtol=ConvX)	# check upper and lower limits !!!
	elif rootf == 'fixed_point':
		if eps == 0.0: n = 0.5
		else: 			
			eqn1 = lambda x: ElectronDensity(params_F,x,mu,Sigman_F,Sigmaa_F,En_F)
			n = fixed_point(eqn1,n_old,xtol=ConvX)
		eqn2 = lambda x: CooperPairDensity(params_F,n,x,Sigman_F,Sigmaa_F,En_F)
		mu = fixed_point(eqn2,mu_old,xtol=ConvX)
		if chat: print '# - - n ={0: .5f}, mu ={1: .5f}'.format(float(n),float(mu))
	if chat: print '# - - n ={0: .5f}, mu ={1: .5f}'.format(float(n),float(mu))
	hfe = ed + U*n


#[GFn_F,GFa_F,EdgePos1,EdgePos2,ABSposGF1,ABSposGF2] = FillGreenHF(U,Delta,GammaR,GammaL,hfe,Phi,mu,wzero,En_F)
#if WriteOut == 1: WriteFile(En_F,GFn_F,GFa_F,params_F,wzero,'HF_green_sc',EminExp,10)

# interacting Green's function ############################
[GFn_F,GFa_F,Det_F,ABS_F,ABSpos_F,Res_F] = FillGreensFunction2(params_F,n,mu,Sigman_F,Sigmaa_F,En_F)
GFn_F = FitTail(En_F,GFn_F,FitMin,FitMax,'odd')
GFa_F = FitTail(En_F,GFa_F,FitMin,FitMax,'even')
wzeroInt = ABS_F[1]
ABSposInt1 = ABSpos_F[0]
GFn_F = FitTail(En_F,GFn_F,FitMin,FitMax,'odd')
GFa_F = FitTail(En_F,GFa_F,FitMin,FitMax,'even')
if WriteOut == 1: WriteFile(En_F,GFn_F,GFa_F,params_F,wzeroInt,'2nd_green',EminExp,10)

###########################################################
# saerching for extremes to see if there is a Hubbard sattelite
'''
from scipy.interpolate import InterpolatedUnivariateSpline
FMin = EdgePos2+1
FMax = FMin + int(15.0/dE)
print '# searching for extremes of GF on ',En_F[FMin],':',En_F[FMax]
Gn  = InterpolatedUnivariateSpline(En_F[FMin:FMax:10],-sp.imag(GFn_F[FMin:FMax:10]))
DGn_F = sp.zeros(len(En_F[FMin:FMax:10]))
for i in range(len(DGn_F)):
	DGn_F[i] = Gn.derivatives(En_F[FMin+10*i])[1]
	#print En_F[FMin+10*i],-sp.imag(GFn_F[FMin+10*i]),Gn(En_F[FMin+10*i]),DGn_F[i]
DGn = InterpolatedUnivariateSpline(En_F[FMin:FMax:10],DGn_F)
Ext_F = DGn.roots()
for i in range(len(Ext_F)):
	if Ext_F[i] > Delta+2.0*(wzero+dE): # minimum due to bubble convolution
		print '# {0: .6f}\t{1: .6f}'.format(float(Ext_F[i]),float(Gn.derivatives(Ext_F[i])[2]))

###########################################################
'''
# final checks and writing the results ####################
n_final  = ElectronDensity(params_F,n,mu,Sigman_F,Sigmaa_F,En_F)
mu_final = CooperPairDensity(params_F,n,mu,Sigman_F,Sigmaa_F,En_F)
if chat: print '# - final densities: n ={0: .5f}, mu ={1: .5f}'.format(float(n_final),float(mu_final))
IDout = IntDOS(GFn_F,En_F)
if chat: print '# - output from 2nd order: ABS = {0: .5f}, intDOS = {1: .5f}'\
.format(float(sp.fabs(wzeroInt)),float(IDout))
if chat: print '# - residues: Gn: [{0: .5f}, {1: .5f}], Ga: {2: .5f}'\
.format(float(Res_F[0]),float(Res_F[1]),float(Res_F[2]))
if sp.sign(Res_F[2]*ResOld_F[2]) < 0.0 and chat:
	print '# Warning: Residue of anomalous function changes sign, {0: .5f} -> {1: .5f}'\
	.format(float(ResOld_F[2]),float(Res_F[2]))
if sp.fabs(Res_F[0]-Res_F[1]) > 1e-3 and eps == 0.0 and chat:
	print "# Warning: Residues of normal GF at ABS don't match."
SEnABS = sp.real(Sigman_F[ABSpos_F[0]])
SEaABS = sp.real(Sigmaa_F[ABSpos_F[0]])
if chat: print '# self-energies at ABS: SEn = {0: .5f}, SEa = {1: .5f}'.format(float(SEnABS),float(SEaABS))

print '{0: .3f}\t{1: .3f}\t{2: .3f}\t{3: .3f}\t{4: .3f}\t{5: .5f}\t{6: .5f}\t{7: .5f}\t{8: .5f}\t{9: .5f}\t{10: .5f}\t{11: .5f}'\
.format(U,GammaR,GammaL,eps,P,float(wzeroInt),float(n),float(mu),Res_F[0],Res_F[1],Res_F[2],float(wzero))

