# SQUAD - superconducting quantum dot
# functions library for case with two sc and one normal lead
# uses scipy, optimized on Python 2.7.5
# Vladislav Pokorny, 2016; pokornyv@fzu.cz

import scipy as sp
from squadlib1 import FillEnergies,SFunctionBand,SFunctionGap,DeltaFunctionBand,DeltaFunctionGap

#####################################################################
# Green's function determinants on real axis ########################

def DetBand(U,GammaR,GammaL,GammaN,Delta,Phi,hfe,mu,x):
	"""	determinant of the HF Green's function in the band region (-inf:-DeltaMax)	"""
	SF = lambda x: SFunctionBand(GammaR,GammaL,Delta,x)
	DF = lambda x: DeltaFunctionBand(GammaR,GammaL,Delta,Phi,x)
	return x**2*(1.0+SF(x))**2-(hfe-1.0j*GammaN/2.0)**2-(DF(x)-U*mu)*(sp.conj(DF(x))-U*mu) #for sp.conj(Delta)!!!
	#return (x*(1.0+SF(x))+GammaN/2.0)**2-hfe**2-(DF(x)-U*mu)**2

def DetGap(U,GammaR,GammaL,GammaN,Delta,Phi,hfe,mu,x):
	"""	determinant of the HF Green's function in the fully gapped region (-DeltaMin:0)	
	    real part only, residues at ABS must be added by hand """
	SF = lambda x: SFunctionGap(GammaR,GammaL,Delta,x)
	DF = lambda x: DeltaFunctionGap(GammaR,GammaL,Delta,Phi,x)
	return x**2*(1.0+SF(x))**2-(hfe-1.0j*GammaN/2.0)**2-(DF(x)-U*mu)*(sp.conj(DF(x))-U*mu) #for sp.conj(Delta)!!!
	#return (x*(1.0+SF(x))+GammaN/2.0)**2-hfe**2-(DF(x)-U*mu)**2 

#####################################################################
# Andreev bound states frequencies ##################################

def AndreevEnergy(U,GammaR,GammaL,GammaN,Delta,Phi,hfe,mu):
	""" returns the ABS frequency in the Hartree-Fock approximation """
	from scipy.optimize import fixed_point	
	SF = lambda x: SFunctionGap(GammaR,GammaL,Delta,x)
	DF = lambda x: DeltaFunctionGap(GammaR,GammaL,Delta,Phi,x)
	eqn = lambda x: sp.sqrt((hfe-1.0j*GammaN/2.0)**2+(DF(x)-U*mu)*(sp.conj(DF(x))-U*mu))/(1.0+SF(x))
	wzero = sp.real_if_close(fixed_point(eqn,0.9*Delta))
	if sp.imag(wzero) != 0.0:
		print "# Warning: non-zero Im w0 = "+str(sp.imag(wzero))
		#wzero = sp.real(wzero)
	return wzero	#fixed_point returns numpy.ndarray!

#####################################################################
# Matsubara sums for HF calculations ################################

def MSum1(U,Delta,GammaR,GammaL,GammaN,hfe,Phi,mu,wzero,En_F):
	"""	Matsubara sum of 1/D(iw) """
	Det = lambda x: DetBand(U,GammaR,GammaL,GammaN,Delta,Phi,hfe,mu,x)
	Det_F = Det(En_F)
	Int_F = sp.imag(Det_F)/(Det_F*sp.conj(Det_F))
	Tail = -Int_F[0]*En_F[0]/2.0	# behaves as -1/x^3 CHECK!!!!
	ContinuumTerm = (simps(Int_F,En_F)+Tail)/sp.pi
	#AndreevTerm = 1.0/DetDiff(U,GammaR,GammaL,GammaN,Delta,Phi,hfe,mu,-wzero)


	return sp.real_if_close(ContinuumTerm + AndreevTerm)

def MSum2(U,Delta,GammaR,GammaL,GammaN,hfe,Phi,mu,wzero,En_F):
	"""	Matsubara sum of (iw(1+s(iw))+eps+GammaN/2)/D(iw) """
	Det = lambda x: DetBand(U,GammaR,GammaL,GammaN,Delta,Phi,hfe,mu,x)
	Sf  = lambda x: SFunctionBand(GammaR,GammaL,Delta,x)
	Det_F = Det(En_F)
	Sf_F =  Sf(En_F)
	Int_F = En_F*(sp.imag(Sf_F)*sp.real(Det_F)-sp.imag(Det_F))/(Det_F*sp.conj(Det_F))
	Tail = -Int_F[0]*En_F[0]	# behaves as 1/x^2 CHECK!!!!
	ContinuumTerm = -(simps(Int_F,En_F)+Tail)/sp.pi
	#AndreevTerm = -wzero*(1.0+SFunctionGap(GammaR,GammaL,GammaN,Delta,-wzero))\
	#/DetDiff(U,GammaR,GammaL,Delta,Phi,hfe,mu,-wzero)


	return sp.real_if_close(ContinuumTerm + AndreevTerm)

def MSum3(U,Delta,GammaR,GammaL,GammaN,hfe,Phi,mu,wzero,En_F):
	"""	Matsubara sum of sp.conj(Delta_Phi(iw))/D(iw) """
	Det = lambda x: DetBand(U,GammaR,GammaL,GammaN,Delta,Phi,hfe,mu,x)
	Del = lambda x: DeltaFunctionBand(GammaR,GammaL,Delta,Phi,x)
	Delta_F = Del(En_F)
	Det_F = Det(En_F)
	#Int_F = -sp.imag(Delta_F)*sp.real(Det_F)/(Det_F*sp.conj(Det_F))	# for sp.conj(Delta)!!!
	Int_F = sp.imag(Delta_F)*sp.real(Det_F)/(Det_F*sp.conj(Det_F))
	Tail = -Int_F[0]*En_F[0]/2.0	# behaves as 1/x^3 CHECK!!!!
	ContinuumTerm = -(simps(Int_F,En_F)+Tail)/sp.pi
	#AndreevTerm = DeltaFunctionGap(GammaR,GammaL,GammaN,Delta,Phi,-wzero)\
	#/DetDiff(U,GammaR,GammaL,Delta,Phi,hfe,mu,-wzero)


	return sp.real_if_close(ContinuumTerm + AndreevTerm)

#####################################################################
# Hartree-Fock solver ###############################################

def SolveHF(U,Delta,GammaR,GammaL,GammaN,eps,P):
	""" HF equations solver """
	Conv = 1e-6			# convergence criterium
	ed = eps-U/2.0		# localized energy level shifted to symmetry point
	ErrMsg = 0			# error message indicator
	imax = 1000			# number of iterations for breaking the cycle
	# filling the arrays ######################################
	dE = 1e-4							# band energy sampling
	Emin = -100.0						# lower cutoff for band energy
	En_F = FillEnergies(Emin,-Delta,dE)	# [Emin:-Delta)
	# initial conditions ######################################
	from scipy.optimize import fixed_point
	Phi = P*sp.pi
	Gamma = (GammaR + GammaL)
	n_density = lambda x: 0.5 - sp.arctan((ed+U*x)/Gamma)/sp.pi	# starting from symmetric metallic case
	n = fixed_point(n_density,0.5)
	mu = 0.2
	hfe = ed+U*n
	wzero = AndreevEnergy(U,GammaR,GammaL,GammaN,Delta,Phi,hfe,mu)
	if Delta == 0:	# HF solution for SIAM ############################
		mu = 0.0
		wzero = 0.0
	else:
		n_old = 1e5
		mu_old = 1e5
		wzero_old = 1e5
		i = 0			# iterations counter
		while any([sp.fabs(sp.real(wzero-wzero_old)) > Conv,sp.fabs(mu-mu_old) > Conv,sp.fabs(n-n_old) > Conv]):
			n_old = n
			mu_old = mu
			wzero_old = wzero
			hfe = ed+U*n
			D1 = MSum1(U,Delta,GammaR,GammaL,GammaN,hfe,Phi,mu,wzero,En_F)
			D3 = MSum3(U,Delta,GammaR,GammaL,GammaN,hfe,Phi,mu,wzero,En_F)
			#print D1,D3
			mu = -D3/(1.0-U*D1)
			if eps == 0: n = 0.5
			else:
				D2 = MSum2(U,Delta,GammaR,GammaL,GammaN,hfe,Phi,mu,wzero,En_F)
				n = (D2+ed*D1)/(1.0-U*D1)
			hfe = ed+U*n
			wzero = AndreevEnergy(U,GammaR,GammaL,GammaN,Delta,Phi,hfe,mu)
			if i > imax: 
				print "# Warning: SolveHF: No convergence after "+str(i)+" iterations, exiting."
				n = mu = wzero = -1.0
				ErrMsg = 1	
				break
			#print '#',i,n,mu,wzero
			i=i+1
	return sp.array([n,mu,wzero,ErrMsg])

