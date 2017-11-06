# SQUAD - superconducting quantum dot
# functions library for general case of asymmetric couplings
# uses scipy, optimized on Python 2.7.5
# Vladislav Pokorny, 2012-2016; pokornyv@fzu.cz

import scipy as sp
from scipy.integrate import simps

#####################################################################
# general functions #################################################

HeavisideTheta = lambda x: (sp.sign(x)+1.0)/2.0
FermiDirac     = lambda E,T: 1.0/(sp.exp((E+1e-12)/T)+1.0)
BoseEinstein   = lambda E,T: 1.0/(sp.exp((E+1e-12)/T)-1.0)


def KondoTemperature(U,Gamma,eps):
	"""	calculating Kondo temperature sqrt(U*Gamma/2)*exp[pi*abs(U**2-4*eps**2)/(8*U*Gamma)] """
	if U != 0.0 and Gamma != 0.0:
		return sp.sqrt(U*Gamma/2.0)*sp.exp(-sp.pi*sp.fabs(U**2-4.0*eps**2)/(8.0*U*Gamma))
	else:
		print "# Warning: KondoTemperature: Kondo temperature not defined."
		return -1.0


def QParticleResidue(En_F,SE_F):
	""" calculates the quasiparticle residue Z = m/m* """
	from scipy.interpolate import InterpolatedUnivariateSpline
	ReSE = InterpolatedUnivariateSpline(En_F,sp.real(SE_F))
	dReSEdw = ReSE.derivatives(0.0)[1]
	return 1.0/(1.0-dReSEdw)


def FillEnergiesLinear(Emin,Emax,dE):
	"""	return the array of energies [Emin,Emin+dE,...,Emax-dE,Emax) """
	En_F = sp.arange(Emin,Emax,dE)
	En_F = sp.around(En_F,int(-sp.log10(dE)))
	return En_F


def FillEnergies(dE,N):
	"""	return the symmetric array of energies 
	[Emin,Emin+dE,...,0,...,-Emin-dE,-Emin] of length N """
	dE_dec = int(-sp.log10(dE))
	En_F = sp.linspace(-(N-1)/2*dE,(N-1)/2*dE,N)
	return sp.around(En_F,dE_dec+2)


def FillEnergies2(dE,N):
	"""	return the symmetric array of energies [Emin,Emin+dE,...,-Emin-dE,-Emin] """
	dE_dec = int(-sp.log10(dE))
	En_F = sp.concatenate([sp.linspace(-N*dE,-dE,N),sp.linspace(0.0,N*dE,N+1)])
	return sp.around(En_F,dE_dec+2)


def FindInEnergies(En_F,w):
	""" return the positions of w in array En_F """
	dE = sp.around(En_F[1]-En_F[0],8)
	dE_dec = int(-sp.log10(dE))
	try:
		Pos = sp.nonzero(En_F == sp.around(w,dE_dec))[0][0]
	except IndexError:
		print 'Warning: w not in En_F'
		Pos = 0
	return Pos

def FindEdges(En_F,Delta):
	""" return the positions of gap edges \pm\Delta in array En_F
	if Delta lies outside, return first and last value """
	dE = sp.around(En_F[1]-En_F[0],8)
	dE_dec = int(-sp.log10(dE))
	if Delta < sp.fabs(En_F[0]):
		EdgePos1 = sp.nonzero(En_F == sp.around(-Delta,dE_dec))[0][0]
		EdgePos2 = sp.nonzero(En_F == sp.around( Delta,dE_dec))[0][0]
	else:
		EdgePos1 = 0
		EdgePos2 = len(En_F)-1
	return sp.array([EdgePos1,EdgePos2])


def Brillouin(J,x):
	"""	Brillouin function """
	a = (2.0*J+1.0)/(2.0*J)
	b = 1.0/(2.0*J)
	return 0.0 if x == 0.0 else a/sp.tanh(a*x)-b/sp.tanh(x)


def Langevin(x):
	""" Langevin function: classical limit of the J=1/2 Brillouin function """
	return 0.0 if x == 0.0 else 1.0/sp.tanh(x)-1.0/x


#####################################################################
# dot-lead hybridizations ###########################################

def SFunctionBand(GammaR,GammaL,Delta,x):
	""" normal hybridization in band region (imaginary) """
	return 1.0j*sp.sign(x)*(GammaL+GammaR)/sp.sqrt(x**2-Delta**2)


def SFunctionGap(GammaR,GammaL,Delta,x):
	""" normal hybridization in gap region (real) """
	return (GammaL+GammaR)/sp.sqrt(Delta**2-x**2)


def DeltaFunctionBand(GammaR,GammaL,Delta,Phi,x):
	""" anomalous hybridization in band region (imaginary)
	PhiC angle helps to keep hybridization real or pure imaginary, not complex
	the 1e-12 offset allows to calculate normal solution w/o superconduvtivity """
	PhiC = sp.arctan((GammaL-GammaR)/(GammaL+GammaR+1e-12)*sp.tan(Phi/2.0))
	return 1.0j*sp.sign(x)*Delta*sp.exp(1.0j*PhiC)/sp.sqrt(x**2-Delta**2)\
	*(GammaL*sp.exp(-1.0j*Phi/2.0) + GammaR*sp.exp(1.0j*Phi/2.0))


def DeltaFunctionGap(GammaR,GammaL,Delta,Phi,x):
	""" anomalous hybridization in gap region (real) """
	PhiC = sp.arctan((GammaL-GammaR)/(GammaL+GammaR+1e-12)*sp.tan(Phi/2.0))
	return Delta*sp.exp(1.0j*PhiC)/sp.sqrt(Delta**2-x**2)\
	*(GammaL*sp.exp(-1.0j*Phi/2.0) + GammaR*sp.exp(1.0j*Phi/2.0))


def SFunctionGapDiff(GammaR,GammaL,Delta,x):
	""" energy derivative of S(w) """
	return x/(Delta**2-x**2)**(3.0/2.0)*(GammaL+GammaR)


def DeltaFunctionGapDiff(GammaR,GammaL,Delta,Phi,x):
	""" energy derivative of Delta(w) """
	PhiC = sp.arctan((GammaL-GammaR)/(GammaL+GammaR+1e-12)*sp.tan(Phi/2.0))	
	return x*Delta*sp.exp(1.0j*PhiC)/(Delta**2-x**2)**(3.0/2.0)\
	*(GammaL*sp.exp(-1.0j*Phi/2.0) + GammaR*sp.exp(1.0j*Phi/2.0))


#####################################################################
# Andreev bound states frequencies ##################################

def AndreevEnergy(U,GammaR,GammaL,Delta,Phi,hfe,mu):
	""" returns the ABS frequency in the Hartree-Fock approximation """
	from scipy.optimize import fixed_point
	SF = lambda x: SFunctionGap(GammaR,GammaL,Delta,x)
	DF = lambda x: DeltaFunctionGap(GammaR,GammaL,Delta,Phi,x)
	eqn = lambda x: sp.sqrt(hfe**2+(DF(x)-U*mu)**2)/(1.0+SF(x))
	## change the initial condition if convercence problems raise
	wzero = sp.real_if_close(fixed_point(eqn,0.999*Delta)) 
	if sp.imag(wzero) != 0.0:
		print "# Warning: non-zero Im w0 = "+str(sp.imag(wzero))
		wzero = sp.real(wzero)
	return wzero	#fixed_point returns numpy.ndarray!


#####################################################################
# Green's function determinants on real axis ########################

def DetBand(U,GammaR,GammaL,Delta,Phi,hfe,mu,x):
	"""	determinant of the HF Green's function in the band region (-inf:-DeltaMax)	"""
	SF = lambda x: SFunctionBand(GammaR,GammaL,Delta,x)
	DF = lambda x: DeltaFunctionBand(GammaR,GammaL,Delta,Phi,x)
	#return x**2*(1.0+SF(x))**2-hfe**2-(DF(x)-U*mu)*(sp.conj(DF(x))-U*mu) #for sp.conj(Delta)!!!
	return x**2*(1.0+SF(x))**2-hfe**2-(DF(x)-U*mu)**2


def DetGap(U,GammaR,GammaL,Delta,Phi,hfe,mu,x):
	"""	determinant of the HF Green's function in the fully gapped region (-DeltaMin:0)	
	    real part only, residues at ABS must be added by hand """
	SF = lambda x: SFunctionGap(GammaR,GammaL,Delta,x)
	DF = lambda x: DeltaFunctionGap(GammaR,GammaL,Delta,Phi,x)
	return x**2*(1.0+SF(x))**2-hfe**2-(DF(x)-U*mu)*(sp.conj(DF(x))-U*mu) #for sp.conj(Delta)!!!
	#return x**2*(1.0+SF(x))**2-hfe**2-(DF(x)-U*mu)**2 


def DetDiff(U,GammaR,GammaL,Delta,Phi,hfe,mu,x):
	""" frequency derivative of the determinant of the HF Green's function in gap region """
	SF = lambda x: SFunctionGap(GammaR,GammaL,Delta,x)
	DF = lambda x: DeltaFunctionGap(GammaR,GammaL,Delta,Phi,x)
	SFD = lambda x: SFunctionGapDiff(GammaR,GammaL,Delta,x)
	DFD = lambda x: DeltaFunctionGapDiff(GammaR,GammaL,Delta,Phi,x)
	#return 2.0*x*(1.0+SF(x))**2+2.0*x**2*(1.0+SF(x))*SFD(x)\
	#-sp.conj(DFD(x))*(DF(x)-U*mu)-DFD(x)*(sp.conj(DF(x))-U*mu) #for sp.conj(Delta)!!!
	DDf = 2.0*x*(1.0+SF(x))**2+2.0*x**2*(1.0+SF(x))*SFD(x)\
	-DFD(x)*(DF(x)-U*mu)-DFD(x)*(DF(x)-U*mu)
	if DDf == 0.0: # some resonance between terms
		print '# Warning: DetDiff: dD/dw = 0, using value 0.01.'
		DDf = 1e-2
	return DDf


def GFnBand(U,GammaR,GammaL,Delta,Phi,hfe,mu,x):
	""" normal GF in band region """
	SF = lambda x: SFunctionBand(GammaR,GammaL,Delta,x)
	DF = lambda x: DeltaFunctionBand(GammaR,GammaL,Delta,Phi,x)
	Det = lambda x: DetBand(U,GammaR,GammaL,Delta,Phi,hfe,mu,x)
	return (x*(1.0+SF(x))+hfe)/Det(x)


def GFnGap(U,GammaR,GammaL,Delta,Phi,hfe,mu,x):
	""" normal GF in gap region """
	SF = lambda x: SFunctionGap(GammaR,GammaL,Delta,x)
	DF = lambda x: DeltaFunctionGap(GammaR,GammaL,Delta,Phi,x)
	Det = lambda x: DetGap(U,GammaR,GammaL,Delta,Phi,hfe,mu,x)
	return (x*(1.0+SF(x))+hfe)/Det(x)


def GFaBand(U,GammaR,GammaL,Delta,Phi,hfe,mu,x):
	""" anomalous GF in band region """
	SF = lambda x: SFunctionBand(GammaR,GammaL,Delta,x)
	DF = lambda x: DeltaFunctionBand(GammaR,GammaL,Delta,Phi,x)
	Det = lambda x: DetBand(U,GammaR,GammaL,Delta,Phi,hfe,mu,x)
	#return -(sp.conj(DF(x))-U*mu)/Det(x) #for sp.conj(Delta)!!!
	return -(DF(x)-U*mu)/Det(x)


def GFaGap(U,GammaR,GammaL,Delta,Phi,hfe,mu,x):
	""" anomalous GF in gap region """
	SF = lambda x: SFunctionGap(GammaR,GammaL,Delta,x)
	DF = lambda x: DeltaFunctionGap(GammaR,GammaL,Delta,Phi,x)
	Det = lambda x: DetGap(U,GammaR,GammaL,Delta,Phi,hfe,mu,x)
	#return -(sp.conj(DF(x))-U*mu)/Det(x) #for sp.conj(Delta)!!!
	return -(DF(x)-U*mu)/Det(x)


def GFresidues(U,Delta,GammaR,GammaL,hfe,Phi,mu,wzero):
	"""	residues of the Green's functions at ABS
	returns an array of three residues: normal electron, normal hole, anomalous	"""
	if Delta == 0.0:
		return sp.array([0.0,0.0,0.0])
	else:
		NomNp =  wzero*(1.0+SFunctionGap(GammaR,GammaL,Delta,wzero))+hfe
		NomNh =  wzero*(1.0+SFunctionGap(GammaR,GammaL,Delta,wzero))-hfe
		NomA  = -(DeltaFunctionGap(GammaR,GammaL,Delta,Phi,wzero)-U*mu)
		Den   =  DetDiff(U,GammaR,GammaL,Delta,Phi,hfe,mu,wzero)
		return sp.real(sp.array([NomNp,NomNh,NomA])/Den)


#####################################################################
# Matsubara sums for HF calculations ################################

def MSum1(U,Delta,GammaR,GammaL,hfe,Phi,mu,wzero,En_F):
	"""	Matsubara sum of 1/D(iw) """
	Det = lambda x: DetBand(U,GammaR,GammaL,Delta,Phi,hfe,mu,x)
	Det_F = Det(En_F)
	Int_F = sp.imag(Det_F)/(Det_F*sp.conj(Det_F))
	Tail = -Int_F[0]*En_F[0]/2.0	# behaves as -1/x^3 CHECK!!!!
	ContinuumTerm = (simps(Int_F,En_F)+Tail)/sp.pi
	AndreevTerm = 1.0/DetDiff(U,GammaR,GammaL,Delta,Phi,hfe,mu,-wzero)
	return sp.real_if_close(ContinuumTerm + AndreevTerm)


def MSum2(U,Delta,GammaR,GammaL,hfe,Phi,mu,wzero,En_F):
	"""	Matsubara sum of (iw(1+s(iw))+eps)/D(iw) """
	Det = lambda x: DetBand(U,GammaR,GammaL,Delta,Phi,hfe,mu,x)
	Sf  = lambda x: SFunctionBand(GammaR,GammaL,Delta,x)
	Det_F = Det(En_F)
	Sf_F =  Sf(En_F)
	Int_F = En_F*(sp.imag(Sf_F)*sp.real(Det_F)-sp.imag(Det_F))/(Det_F*sp.conj(Det_F))
	Tail = -Int_F[0]*En_F[0]	# behaves as 1/x^2 CHECK!!!!
	ContinuumTerm = -(simps(Int_F,En_F)+Tail)/sp.pi
	AndreevTerm = -wzero*(1.0+SFunctionGap(GammaR,GammaL,Delta,-wzero))\
	/DetDiff(U,GammaR,GammaL,Delta,Phi,hfe,mu,-wzero)
	return sp.real_if_close(ContinuumTerm + AndreevTerm)


def MSum3(U,Delta,GammaR,GammaL,hfe,Phi,mu,wzero,En_F):
	"""	Matsubara sum of sp.conj(Delta_Phi(iw))/D(iw) """
	Det = lambda x: DetBand(U,GammaR,GammaL,Delta,Phi,hfe,mu,x)
	Del = lambda x: DeltaFunctionBand(GammaR,GammaL,Delta,Phi,x)
	Delta_F = Del(En_F)
	Det_F = Det(En_F)
	#Int_F = -sp.imag(Delta_F)*sp.real(Det_F)/(Det_F*sp.conj(Det_F))	# for sp.conj(Delta)!!!
	Int_F = sp.imag(Delta_F)*sp.real(Det_F)/(Det_F*sp.conj(Det_F))
	Tail = -Int_F[0]*En_F[0]/2.0	# behaves as 1/x^3 CHECK!!!!
	ContinuumTerm = -(simps(Int_F,En_F)+Tail)/sp.pi
	AndreevTerm = DeltaFunctionGap(GammaR,GammaL,Delta,Phi,-wzero)\
	/DetDiff(U,GammaR,GammaL,Delta,Phi,hfe,mu,-wzero)
	return sp.real_if_close(ContinuumTerm + AndreevTerm)


def MSumsHF(U,Delta,GammaR,GammaL,hfe,Phi,mu,wzero,En_F):
	"""	Matsubara sum of 1/D(iw)
		Matsubara sum of (iw(1+s(iw))+eps)/D(iw)
		Matsubara sum of sp.conj(Delta_Phi(iw))/D(iw) """
	Det = lambda x: DetBand(U,GammaR,GammaL,Delta,Phi,hfe,mu,x)
	Sf  = lambda x: SFunctionBand(GammaR,GammaL,Delta,x)
	Del = lambda x: DeltaFunctionBand(GammaR,GammaL,Delta,Phi,x)
	Det_F = Det(En_F)
	Sf_F =  Sf(En_F)
	Delta_F = Del(En_F)
	if len(En_F)>0:	# Delta smaller than energy minimum
		Int1_F = sp.imag(Det_F)/(Det_F*sp.conj(Det_F))
		Int2_F = En_F*(sp.imag(Sf_F)*sp.real(Det_F)-sp.imag(Det_F))/(Det_F*sp.conj(Det_F))
		Int3_F = sp.imag(Delta_F)*sp.real(Det_F)/(Det_F*sp.conj(Det_F))
		Tail1 = -Int1_F[0]*En_F[0]/2.0	# behaves as -1/x^3 CHECK!!!!
		Tail2 = -Int2_F[0]*En_F[0]	    # behaves as 1/x^2 CHECK!!!!
		Tail3 = -Int3_F[0]*En_F[0]/2.0	# behaves as 1/x^3 CHECK!!!!
		ContTerm1 = (simps(Int1_F,En_F)+Tail1)/sp.pi
		ContTerm2 = -(simps(Int2_F,En_F)+Tail2)/sp.pi
		ContTerm3 = -(simps(Int3_F,En_F)+Tail3)/sp.pi
	else:
		ContTerm1 = ContTerm2 = ContTerm3 = 0.0
	AndreevTerm1 = 1.0/DetDiff(U,GammaR,GammaL,Delta,Phi,hfe,mu,-wzero)
	AndreevTerm2 = -wzero*(1.0+SFunctionGap(GammaR,GammaL,Delta,-wzero))/DetDiff(U,GammaR,GammaL,Delta,Phi,hfe,mu,-wzero)
	AndreevTerm3 = DeltaFunctionGap(GammaR,GammaL,Delta,Phi,-wzero)/DetDiff(U,GammaR,GammaL,Delta,Phi,hfe,mu,-wzero)
	D1 = sp.real_if_close(ContTerm1 + AndreevTerm1)
	D2 = sp.real_if_close(ContTerm2 + AndreevTerm2)
	D3 = sp.real_if_close(ContTerm3 + AndreevTerm3)
	return [D1,D2,D3]


#####################################################################
# Hartree-Fock solver ###############################################

def SolveHF(Params_F):
	""" HF equations solver """
	from scipy.optimize import fixed_point
	import params_ss as pss
	[U,Delta,GammaR,GammaL,GammaN,Phi,eps,h] = Params_F
	Gamma = (GammaR + GammaL)
	ed = eps-U/2.0		# localized energy level shifted to symmetry point
	ErrMsg = 0			# error message indicator
	## filling the arrays #################################
	dE = 1e-4							# band energy sampling
	Emin = -100.0						# lower cutoff for band energy
	En_F = FillEnergiesLinear(Emin,-Delta,dE)	# [Emin:-Delta)
	## initial conditions #################################
	## change these if no convergence is achieved #########
	n_density = lambda x: 0.5 - sp.arctan((ed+U*x)/Gamma)/sp.pi	# starting from symmetric metallic case
	n = fixed_point(n_density,0.5)
	#n = 0.001
	mu = 0.01
	hfe = ed+U*n
	wzero = AndreevEnergy(U,GammaR,GammaL,Delta,Phi,hfe,mu)
	if Delta == 0:	# HF solution for SIAM ############################
		mu = 0.2
		wzero = 0.0
	else:
		n_old = 1e5
		mu_old = 1e5
		wzero_old = 1e5
		i = 0			            # iterations counter
		imax = pss.P['HF_max_iter']	# number of iterations for breaking the cycle
		while any([sp.fabs(wzero-wzero_old) > pss.P['ConvHF'],\
		           sp.fabs(mu-mu_old)       > pss.P['ConvHF'],\
                   sp.fabs(n-n_old)         > pss.P['ConvHF']]):
			n_old = n
			mu_old = mu
			wzero_old = wzero
			hfe = ed+U*n
			[D1,D2,D3] = MSumsHF(U,Delta,GammaR,GammaL,hfe,Phi,mu,wzero,En_F)
			mu = -D3/(1.0-U*D1)
			if eps == 0: n = 0.5
			else:
				n = (D2+ed*D1)/(1.0-U*D1)
			hfe = ed+U*n
			wzero = AndreevEnergy(U,GammaR,GammaL,Delta,Phi,hfe,mu)
			#print n,mu,wzero
			if i > imax: 
				print "# Warning: SolveHF: No convergence after "+str(i)+" iterations, exit."
				n = mu = wzero = -1.0
				ErrMsg = 1	
				break
			if sp.imag(mu) != 0.0:
				print "# Warning: non-zero Im mu = "+str(sp.imag(mu))
				mu = sp.real(mu)
			i=i+1
	if pss.P['WriteIO']: print('#   {0: 3d} iterations,\t n = {1: .6f} +{2: .6f}i,  mu = {3: .6f} +{4: .6f}i'\
	.format(i,float(sp.real(n)),float(sp.imag(n)),float(sp.real(mu)),float(sp.imag(mu))))
	return sp.array([n,mu,wzero,ErrMsg])


#####################################################################
# filling the arrays with HF GFs ####################################

def FillGreenHF(U,Delta,GammaR,GammaL,hfe,Phi,mu,wzero,En_F):
	""" filling the En_F array with HF Green's functions """
	## define the lambdas
	GFn_band = lambda x: GFnBand(U,GammaR,GammaL,Delta,Phi,hfe,mu,x)
	GFa_band = lambda x: GFaBand(U,GammaR,GammaL,Delta,Phi,hfe,mu,x)
	GFn_gap  = lambda x: GFnGap(U,GammaR,GammaL,Delta,Phi,hfe,mu,x)
	GFa_gap  = lambda x: GFaGap(U,GammaR,GammaL,Delta,Phi,hfe,mu,x)
	## find the frequency mesh
	dE = sp.around(En_F[1]-En_F[0],8)
	dE_dec = int(-sp.log10(dE))
	zero_F = sp.array([0.0])	# gap edge
	## find special points
	[EdgePos1,EdgePos2] = FindEdges(En_F,Delta)
	if sp.fabs(wzero) > dE:
		ABSpos1 = sp.nonzero(En_F==sp.around(-wzero,dE_dec))[0][0]
		ABSpos2 = sp.nonzero(En_F==sp.around(wzero,dE_dec))[0][0]
	else:	# putting poles at lowest possible points
		ABSpos1 = sp.nonzero(En_F==sp.around(-dE,dE_dec))[0][0]
		ABSpos2 = sp.nonzero(En_F==sp.around(dE,dE_dec))[0][0]
	## fill the arrays
	GFn_F = sp.concatenate((GFn_band(En_F[:EdgePos1]),zero_F,GFn_gap(En_F[EdgePos1+1:EdgePos2])\
	,zero_F,GFn_band(En_F[EdgePos2+1:])))
	GFa_F = sp.concatenate((GFa_band(En_F[:EdgePos1]),zero_F,GFa_gap(En_F[EdgePos1+1:EdgePos2])\
	,zero_F,GFa_band(En_F[EdgePos2+1:])))
	## calculate residues at ABS
	[ResNp1,ResNh1,ResA1] = GFresidues(U,Delta,GammaR,GammaL,hfe,Phi,mu,-wzero)
	[ResNp2,ResNh2,ResA2] = GFresidues(U,Delta,GammaR,GammaL,hfe,Phi,mu, wzero)
	## add singular parts from ABS to imag. parts
	GFn_F[ABSpos1]=-1.0j*ResNp1*sp.pi/dE
	GFa_F[ABSpos1]=-1.0j*ResA1*sp.pi/dE
	GFn_F[ABSpos2]=-1.0j*ResNp2*sp.pi/dE
	GFa_F[ABSpos2]=-1.0j*ResA2*sp.pi/dE
	return [GFn_F,GFa_F,EdgePos1,EdgePos2,ABSpos1,ABSpos2]

## squadlib1.py end ##

