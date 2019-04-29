################################################################
# SQUAD - superconducting quantum dot                          #
# Copyright (C) 2012-2019  Vladislav Pokorny; pokornyv@fzu.cz  #
# homepage: github.com/pokornyv/SQUAD                          #
# squadlib1.py - library of functions                          #
# method described in                                          #
#    Sci. Rep. 5, 8821 (2015).                                 #
#    Phys. Rev. B 93, 024523 (2016).                           #
################################################################

from config_squad import *
from scipy.integrate import trapz,simps
from scipy.optimize import fixed_point,brentq

#############################
##### List of functions: ####
# KondoTemperature
# FindInEnergies
# SFunctionBand
# SFunctionGap
# DeltaFunctionBand
# DeltaFunctionGap
# SFunctionGapDiff
# DeltaFunctionGapDiff
# AndreevEnergy
# DetBand
# DetGap
# DetDiff
# GFnBand
# GFnGap
# GFaBand
# GFaGap
# GFresidues
# FillGreenHF
# MSumsHF
# SolveHF

#####################################################################
# general functions #################################################

def KondoTemperature():
	""" Kondo temperature from Bethe ansatz for single impurity Anderson model """
	if U == 0.0 or GammaTot == 0.0:
		print('# Warning: KondoTemperature: not defined.')
		TK = -1.0
	else:
		TK = sp.sqrt(U*GammaTot/2.0)*sp.exp(-sp.pi*sp.fabs(U**2-4.0*eps**2)/(8.0*U*GammaTot))
	return TK


def FindInEnergies(x,X_A):
	""" returns positions of points in energy array X_A """
	dX = sp.around(X_A[1]-X_A[0],8)
	dX_dec = int(-sp.log10(dX))
	Pos = sp.nonzero(sp.around(X_A,dX_dec) == sp.around(x,dX_dec))[0][0]
	return Pos

#####################################################################
# dot-lead hybridizations ###########################################

def SFunctionBand(x):
	""" normal hybridization in band region (imaginary) """
	return 1.0j*sp.sign(x)*(GammaL+GammaR)/sp.sqrt(x**2-Delta**2)


def SFunctionGap(x):
	""" normal hybridization in gap region (real) """
	return (GammaL+GammaR)/sp.sqrt(Delta**2-x**2)


def DeltaFunctionBand(x):
	""" anomalous hybridization in band region (imaginary)
	    PhiC angle helps to keep hybridization real or pure imaginary, not complex """
	PhiC = sp.arctan((GammaL-GammaR)/(GammaL+GammaR+1e-12)*sp.tan(Phi/2.0))
	return 1.0j*sp.sign(x)*Delta*sp.exp(1.0j*PhiC)/sp.sqrt(x**2-Delta**2)\
	*(GammaL*sp.exp(-1.0j*Phi/2.0) + GammaR*sp.exp(1.0j*Phi/2.0))


def DeltaFunctionGap(x):
	""" anomalous hybridization in gap region (real) """
	PhiC = sp.arctan((GammaL-GammaR)/(GammaL+GammaR+1e-12)*sp.tan(Phi/2.0))
	return Delta*sp.exp(1.0j*PhiC)/sp.sqrt(Delta**2-x**2)\
	*(GammaL*sp.exp(-1.0j*Phi/2.0) + GammaR*sp.exp(1.0j*Phi/2.0))


def SFunctionGapDiff(x):
	""" energy derivative of S(w) """
	return x/(Delta**2-x**2)**(3.0/2.0)*(GammaL+GammaR)


def DeltaFunctionGapDiff(x):
	""" energy derivative of Delta(w) """
	PhiC = sp.arctan((GammaL-GammaR)/(GammaL+GammaR+1e-12)*sp.tan(Phi/2.0))	
	return x*Delta*sp.exp(1.0j*PhiC)/(Delta**2-x**2)**(3.0/2.0)\
	*(GammaL*sp.exp(-1.0j*Phi/2.0) + GammaR*sp.exp(1.0j*Phi/2.0))


## hybridizations as lambda functions
SFb = lambda x: SFunctionBand(x)
DFb = lambda x: DeltaFunctionBand(x)
SFg = lambda x: SFunctionGap(x)
DFg = lambda x: DeltaFunctionGap(x)
SFD = lambda x: SFunctionGapDiff(x)
DFD = lambda x: DeltaFunctionGapDiff(x)

#####################################################################
# Andreev bound states frequencies ##################################

def AndreevEnergy(hfe,mu):
	""" returns the ABS frequency in the Hartree-Fock approximation """
	eqn = lambda x: sp.sqrt(hfe**2+(DFg(x)-U*mu)**2)/(1.0+SFg(x))
	## change the initial condition init_val if convercence problems raise
	wzero = sp.real_if_close(fixed_point(eqn,ABSinit_val*Delta))
	if sp.fabs(sp.imag(wzero)) > 1e-12:
		print("# - Warning: AndreevEnergy: Non-zero Im w(ABS) = {0: .5e}".format(sp.imag(wzero)))
	return sp.float64(sp.real(wzero))	## caution: fixed_point returns numpy.ndarray

#####################################################################
# Green function determinants #######################################

def DetBand(hfe,mu,x):
	""" determinant of the HF Green function in the band region (-inf:-DeltaMax)	"""
	return x**2*(1.0+SFb(x))**2-hfe**2-(DFb(x)-U*mu)**2


def DetGap(hfe,mu,x):
	""" determinant of the HF Green function in the fully gapped region (-DeltaMin:0)	
	    real part only, residues at ABS must be added by hand """
	return x**2*(1.0+SFg(x))**2-hfe**2-(DFg(x)-U*mu)*(sp.conj(DFg(x))-U*mu) #for sp.conj(Delta)!!!


def DetDiff(hfe,mu,x):
	""" frequency derivative of the determinant of the HF Green function in gap region """
	DDf = 2.0*x*(1.0+SFg(x))**2+2.0*x**2*(1.0+SFg(x))*SFD(x)\
	-DFD(x)*(DFg(x)-U*mu)-DFD(x)*(DFg(x)-U*mu)
	if DDf == 0.0: ## some resonance between terms?
		print('# - Warning: DetDiff: dD/dw = 0, using value 0.01.')
		DDf = 1e-2
	return DDf

#####################################################################
# Green functions ###################################################

def GFnBand(hfe,mu,x):
	""" normal GF in band region """
	return (x*(1.0+SFb(x))+hfe)/DetBand(hfe,mu,x)


def GFnGap(hfe,mu,x):
	""" normal GF in gap region """
	return (x*(1.0+SFg(x))+hfe)/DetGap(hfe,mu,x)


def GFaBand(hfe,mu,x):
	""" anomalous GF in band region """
	return -(DFb(x)-U*mu)/DetBand(hfe,mu,x)


def GFaGap(hfe,mu,x):
	""" anomalous GF in gap region """
	return -(DFg(x)-U*mu)/DetGap(hfe,mu,x)


def GFresidues(hfe,mu,wzero):
	""" residues of the Green functions at ABS
	    returns an array of three residues: normal electron, normal hole, anomalous """
	NomNp =  wzero*(1.0+SFg(wzero))+hfe
	NomNh =  wzero*(1.0+SFg(wzero))-hfe
	NomA  = -(DFg(wzero)-U*mu)
	return sp.real(sp.array([NomNp,NomNh,NomA])/DetDiff(hfe,mu,wzero))


def FillGreenHF(hfe,mu,wzero):
	""" filling the arrays with HF Green functions """
	## define the lambdas for given hfe and mu
	GFn_band = lambda x: GFnBand(hfe,mu,x)
	GFa_band = lambda x: GFaBand(hfe,mu,x)
	GFn_gap  = lambda x: GFnGap(hfe,mu,x)
	GFa_gap  = lambda x: GFaGap(hfe,mu,x)
	## gap edge
	zero_A = sp.array([0.0])
	## find special points
	if sp.fabs(wzero) > dE:
		[ABSpos1,ABSpos2] = [FindInEnergies(-wzero,En_A),FindInEnergies( wzero,En_A)]
	else:	# putting poles at lowest possible points
		print('# - Warning: FillGreenHF: ABS very close to Fermi energy.')
		[ABSpos1,ABSpos2] = [FindInEnergies(-dE,En_A),FindInEnergies(dE,En_A)]
	## fill the arrays
	GFn_A = sp.concatenate((GFn_band(En_A[:EdgePos1]),sp.zeros(1)\
	,GFn_gap(En_A[EdgePos1+1:EdgePos2])\
	,sp.zeros(1),GFn_band(En_A[EdgePos2+1:])))
	GFa_A = sp.concatenate((GFa_band(En_A[:EdgePos1]),sp.zeros(1),\
	GFa_gap(En_A[EdgePos1+1:EdgePos2])\
	,sp.zeros(1),GFa_band(En_A[EdgePos2+1:])))
	## calculate residues at ABS
	[ResNp1,ResNh1,ResA1] = GFresidues(hfe,mu,-wzero)
	[ResNp2,ResNh2,ResA2] = GFresidues(hfe,mu, wzero)
	## add ABS to imag. parts
	GFn_A[ABSpos1]=-1.0j*ResNp1*sp.pi/dE
	GFa_A[ABSpos1]=-1.0j*ResA1 *sp.pi/dE
	GFn_A[ABSpos2]=-1.0j*ResNp2*sp.pi/dE
	GFa_A[ABSpos2]=-1.0j*ResA2 *sp.pi/dE
	return [GFn_A,GFa_A,ABSpos1,ABSpos2]

#####################################################################
# Matsubara sums for HF calculations ################################

def MSumsHF(hfe,mu,wzero,X_A):
	"""	Matsubara sum of 1/Det(iw)
		Matsubara sum of (iw(1+s(iw))+eps)/Det(iw)
		Matsubara sum of sp.conj(Delta(iw))/Det(iw) """
	Det_A   = DetBand(hfe,mu,X_A)
	Int1_A = sp.imag(Det_A)/(Det_A*sp.conj(Det_A))
	Int2_A = X_A*(sp.imag(SFb(X_A))*sp.real(Det_A)-sp.imag(Det_A))/(Det_A*sp.conj(Det_A))
	Int3_A = sp.imag(DFb(X_A))*sp.real(Det_A)/(Det_A*sp.conj(Det_A))
	Tail1 = -Int1_A[0]*X_A[0]/2.0	## behaves as -1/x^3
	Tail2 = -Int2_A[0]*X_A[0]	## behaves as  1/x^2
	Tail3 = -Int3_A[0]*X_A[0]/2.0	## behaves as  1/x^3
	ContTerm1 =  (simps(Int1_A,X_A)+Tail1)/sp.pi
	ContTerm2 = -(simps(Int2_A,X_A)+Tail2)/sp.pi
	ContTerm3 = -(simps(Int3_A,X_A)+Tail3)/sp.pi
	AndreevTerm1 =  1.0/DetDiff(hfe,mu,-wzero)
	AndreevTerm2 = -wzero*(1.0+SFunctionGap(-wzero))/DetDiff(hfe,mu,-wzero)
	AndreevTerm3 =  DeltaFunctionGap(-wzero)/DetDiff(hfe,mu,-wzero)
	D1 = sp.real_if_close(ContTerm1 + AndreevTerm1)
	D2 = sp.real_if_close(ContTerm2 + AndreevTerm2)
	D3 = sp.real_if_close(ContTerm3 + AndreevTerm3)
	return [D1,D2,D3]

#####################################################################
# The Hartree-Fock solver ###########################################

def SolveHF():
	""" Hartree-Fock equations solver """
	ed = eps-U/2.0            ## local energy level shifted to symmetry point
	ErrMsg = 0                ## error message indicator
	## filling the arrays #################################
	dX   = 1e-4         ## band energy sampling
	Xmin = -100.0       ## lower cutoff for band energy
	X_A = sp.arange(Xmin,-Delta,dX)
	## initial conditions #################################
	## change these if no convergence is achieved #########
	n_siam = lambda x: 0.5 - sp.arctan((ed+U*x)/(GammaR+GammaL))/sp.pi	
	n = fixed_point(n_siam,0.5)
	mu = 0.2
	hfe = ed+U*n
	wzero = AndreevEnergy(hfe,mu)
	n_old = 1e5
	mu_old = 1e5
	wzero_old = 1e5
	k = 0
	while any([sp.fabs(wzero-wzero_old) > ConvHF,\
	           sp.fabs(mu-mu_old)       > ConvHF,\
                sp.fabs(n-n_old)         > ConvHF]):
		[n_old,mu_old,wzero_old] = [n,mu,wzero]
		hfe = ed+U*n
		[D1,D2,D3] = MSumsHF(hfe,mu,wzero,X_A)
		mu = -D3/(1.0-U*D1)
		n = 0.5 if eps == 0.0 else (D2+ed*D1)/(1.0-U*D1)
		hfe = ed+U*n
		wzero = AndreevEnergy(hfe,mu)
		if k > HF_max_iter: 
			print('# - Error: SolveHF: No convergence after {0: 5d} iterations, exit.'.format(k))
			n = mu = wzero = -1.0
			ErrMsg = 1
			break
		if sp.fabs(sp.imag(n))  > 1e-12:
			print('# - Warning: SolveHF: neglecting non-zero Im n  = {0: .5e}'.format(sp.imag(n)))
		n = sp.real(n)
		if sp.fabs(sp.imag(mu)) > 1e-12:
			print('# - Warning: SolveHF: neglecting non-zero Im mu = {0: .5e}'.format(sp.imag(mu)))
		mu = sp.real(mu)
		k += 1
	if chat: print('# - Converged after {0: 3d} iterations,  n = {1: .6f},  mu = {2: .6f}'\
	.format(k,float(n),float(mu)))
	return sp.array([n,mu,wzero,ErrMsg])

## squadlib1.py end ##

