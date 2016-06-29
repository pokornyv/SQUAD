# SQUAD - superconducting quantum dot
# functions library for case with two sc and one normal lead
# uses scipy, optimized on Python 2.7.5
# Vladislav Pokorny, 2016; pokornyv@fzu.cz

import scipy as sp
from scipy.integrate import simps
from scipy.fftpack import fft,ifft
from squadlib1 import SFunctionBand,SFunctionGap,DeltaFunctionBand,DeltaFunctionGap
from squadlib1 import FermiDirac,BoseEinstein
from squadlib2 import KramersKronigFFT
import params as p

#####################################################################
# Green's function determinants on real axis ########################

def DetBand(params_F,hfe,mu,x):
	"""	determinant of the HF Green's function in 
	the band region (-inf:-Delta),(Delta,inf) """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps] = params_F 
	SF = lambda x: SFunctionBand(GammaR,GammaL,Delta,x+p.offset_x)
	DF = lambda x: DeltaFunctionBand(GammaR,GammaL,Delta,Phi,x+p.offset_x)
	return (x*(1.0+SF(x))+1.0j*GammaNbar)**2-hfe**2-(DF(x)-U*mu)**2


def DetGap(params_F,hfe,mu,x):
	"""	determinant of the HF Green's function in 
	the gap region (-Delta:Delta) """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps] = params_F 
	SF = lambda x: SFunctionGap(GammaR,GammaL,Delta,x+p.offset_x)
	DF = lambda x: DeltaFunctionGap(GammaR,GammaL,Delta,Phi,x+p.offset_x)
	return (x*(1.0+SF(x))+1.0j*GammaNbar)**2-hfe**2-(DF(x)-U*mu)**2 


def GFnBand(params_F,hfe,mu,x):
	""" normal GF in band region """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps] = params_F 
	SF = lambda x: SFunctionBand(GammaR,GammaL,Delta,x+p.offset_x)
	DF = lambda x: DeltaFunctionBand(GammaR,GammaL,Delta,Phi,x+p.offset_x)
	Det = lambda x: DetBand(params_F,hfe,mu,x)
	return (x*(1.0+SF(x))+1.0j*GammaNbar+hfe)/Det(x)


def GFnGap(params_F,hfe,mu,x):
	""" normal GF in gap region """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps] = params_F 
	SF = lambda x: SFunctionGap(GammaR,GammaL,Delta,x+p.offset_x)
	DF = lambda x: DeltaFunctionGap(GammaR,GammaL,Delta,Phi,x+p.offset_x)
	Det = lambda x: DetGap(params_F,hfe,mu,x)
	return (x*(1.0+SF(x))+1.0j*GammaNbar+hfe)/Det(x)


def GFaBand(params_F,hfe,mu,x):
	""" anomalous GF in band region """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps] = params_F 
	SF = lambda x: SFunctionBand(GammaR,GammaL,Delta,x+p.offset_x)
	DF = lambda x: DeltaFunctionBand(GammaR,GammaL,Delta,Phi,x+p.offset_x)
	Det = lambda x: DetBand(params_F,hfe,mu,x)
	return -(DF(x)-U*mu)/Det(x)


def GFaGap(params_F,hfe,mu,x):
	""" anomalous GF in gap region """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps] = params_F 
	SF = lambda x: SFunctionGap(GammaR,GammaL,Delta,x+p.offset_x)
	DF = lambda x: DeltaFunctionGap(GammaR,GammaL,Delta,Phi,x+p.offset_x)
	Det = lambda x: DetGap(params_F,hfe,mu,x)
	return -(DF(x)-U*mu)/Det(x)


def GreensFunction(params_F,hfe,mu,SEn_F,SEa_F,SEnStar_F,SEaStar_F,En_F,A):
	"""	general form the Nambu Green's function and its determinant
	A='band'/'gap' calculates determinant for band/gap """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps] = params_F 
	#hfe = eps+U*(n-0.5)
	if A == 'band':
		S_F = SFunctionBand(GammaR,GammaL,Delta,En_F+p.offset_x)
		D_F = DeltaFunctionBand(GammaR,GammaL,Delta,Phi,En_F+p.offset_x)
	else:
		S_F = SFunctionGap(GammaR,GammaL,Delta,En_F+p.offset_x)
		D_F = DeltaFunctionGap(GammaR,GammaL,Delta,Phi,En_F+p.offset_x)
	Detn_F = (En_F*(1.0+S_F)+1.0j*GammaNbar-hfe-SEn_F)*(En_F*(1.0+S_F)+1.0j*GammaNbar+hfe-SEnStar_F)
	Deta_F = (D_F-U*mu-SEa_F)*(D_F-U*mu-SEaStar_F)
	Det_F = Detn_F-Deta_F
	GFn_F =  (En_F*(1.0+S_F)+1.0j*GammaNbar+hfe-SEnStar_F)/Det_F
	GFa_F = -(D_F-U*mu-SEa_F)/Det_F
	return [GFn_F,GFa_F,Det_F]


#####################################################################
# Hartree-Fock solver ###############################################

def SolveHFssn(params_F,chat):
	""" HF equations solver for case with two SC and one normal lead """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps] = params_F 
	ed = eps-U/2.0		# localized energy level shifted to symmetry point
	ErrMsg = 0			# error message indicator
	## fill the energy array
	dE = 1e-4							# band energy sampling
	Emin = -100.0						# lower cutoff for band energy
	En_F = sp.around(sp.arange(Emin,0.0,dE),int(-sp.log10(dE)))
	## initial condition: normal state w/o superconductivity
	from scipy.optimize import fixed_point
	Gamma = (GammaR + GammaL + GammaNbar)
	n_density = lambda x: 0.5 - sp.arctan((ed+U*x)/Gamma)/sp.pi
	n = fixed_point(n_density,0.5)
	mu = 0.2	# good guess in zero-phase, pi-phase needs small and negative
	hfe = ed+U*n
	if Delta == 0: mu = 0.0
	else:
		n_old = 1e5
		mu_old = 1e5
		i = 0		# iterations counter
		imax = p.HF_max_iter	# maximum number of iterations
		while any([sp.fabs(mu-mu_old) > p.ConvHF, sp.fabs(n-n_old) > p.ConvHF]):
			n_old = n
			mu_old = mu
			hfe = ed+U*n
			[D1,D2,D3] = MSums(params_F,hfe,mu,En_F)
			#print(D1,D2,D3)
			mu = -D3/(1.0-U*D1)
			if eps == 0: n = 0.5
			else: n = (D2+ed*D1)/(1.0-U*D1)
			hfe = ed+U*n
			if i > imax: 
				print "# Warning: SolveHF: No convergence after "+str(i)+" iterations, exit."
				n = mu = -1.0
				ErrMsg = 1	
				break
			if chat: print('# {0: 3d}\t n = {1: .6f} +{2: .6f}i,  mu = {3: .6f} +{4: .6f}i'\
			.format(i,float(sp.real(n)),float(sp.imag(n)),float(sp.real(mu)),float(sp.imag(mu))))
			i += 1
	return sp.array([n,mu,ErrMsg])


#####################################################################
# convolution procedures using FFT ##################################

def TwoParticleBubbles(GFn_F,GFa_F,En_F):
	""" calculating the two-particle bubbles using fft from scipy.fftpack """
	N = (len(En_F)-1)/2
	dE = sp.around(En_F[1]-En_F[0],8)
	GFnh_F = -sp.conj(sp.flipud(GFn_F))	
	## zero-padding of the arrays to double the size
	FD_F = 1.0*sp.concatenate([sp.zeros(3*N+3),sp.ones(N+1)])
	ImGFp_F = sp.concatenate([sp.imag(GFn_F[N:]),sp.zeros(2*N+3),sp.imag(GFn_F[:N])])
	ImGFh_F = sp.concatenate([sp.imag(GFnh_F[N:]),sp.zeros(2*N+3),sp.imag(GFnh_F[:N])])
	ImGFa_F = sp.concatenate([sp.imag(GFa_F[N:]),sp.zeros(2*N+3),sp.imag(GFa_F[:N])])
	## perform convolution/cross-correlation via FFT 
	ftImChin1_F = -sp.conj(fft(FD_F*ImGFp_F))*fft(ImGFp_F)*dE
	ftImChin2_F = fft(FD_F*ImGFp_F)*sp.conj(fft(ImGFp_F))*dE
	ImChin_F = sp.real(ifft(ftImChin1_F+ftImChin2_F))/sp.pi
	ImChin_F = sp.concatenate([ImChin_F[3*N+4:],ImChin_F[:N+1]])
	ftImChia1_F = -sp.conj(fft(FD_F*ImGFa_F))*fft(ImGFa_F)*dE
	ftImChia2_F = fft(FD_F*ImGFa_F)*sp.conj(fft(ImGFa_F))*dE
	ImChia_F = sp.real(ifft(ftImChia1_F+ftImChia2_F))/sp.pi
	ImChia_F = sp.concatenate([ImChia_F[3*N+4:],ImChia_F[:N+1]])
	## find real part from imaginary using KK relations
	Chin_F = KramersKronigFFT(ImChin_F) + 1.0j*ImChin_F
	Chia_F = KramersKronigFFT(ImChia_F) + 1.0j*ImChia_F
	return [Chin_F,Chia_F]


def SelfEnergy(GFn_F,GFa_F,ChiGamma_F,En_F):
	""" calculating the dynamical self-energy from Schwinger-Dyson equation
	uses fft from scipy.fftpack """
	N = (len(En_F)-1)/2
	dE = sp.around(En_F[1]-En_F[0],8)
	## zero-padding of the arrays to double the size
	FD_F = 1.0*sp.concatenate([sp.zeros(3*N+3),sp.ones(N+1)])
	BE_F = -sp.copy(FD_F)
	ImGFn_F = sp.concatenate([sp.imag(GFn_F[N:]),sp.zeros(2*N+3),sp.imag(GFn_F[:N])])
	ImGFa_F = sp.concatenate([sp.imag(GFa_F[N:]),sp.zeros(2*N+3),sp.imag(GFa_F[:N])])
	ImCG_F = sp.concatenate([sp.imag(ChiGamma_F[N:]),sp.zeros(2*N+3),sp.imag(ChiGamma_F[:N])])
	## perform convolution/cross-correlation via FFT 
	ftImSEn1_F = -sp.conj(fft(BE_F*ImCG_F))*fft(ImGFn_F)*dE
	ftImSEn2_F = -fft(FD_F*ImGFn_F)*sp.conj(fft(ImCG_F))*dE
	ImSEn_F = sp.real(ifft(ftImSEn1_F+ftImSEn2_F))/sp.pi
	ImSEn_F = sp.concatenate([ImSEn_F[3*N+4:],ImSEn_F[:N+1]])
	ftImSEa1_F = -sp.conj(fft(BE_F*ImCG_F))*fft(ImGFa_F)*dE
	ftImSEa2_F = -fft(FD_F*ImGFa_F)*sp.conj(fft(ImCG_F))*dE
	ImSEa_F = sp.real(ifft(ftImSEa1_F+ftImSEa2_F))/sp.pi
	ImSEa_F = sp.concatenate([ImSEa_F[3*N+4:],ImSEa_F[:N+1]])
	## find real part from imaginary using KK relations
	Sigman_F = KramersKronigFFT(ImSEn_F) + 1.0j*ImSEn_F
	Sigmaa_F = KramersKronigFFT(ImSEa_F) + 1.0j*ImSEa_F
	return [Sigman_F,Sigmaa_F]


#####################################################################
# functions to calculate particle densities #########################

def MSums(params_F,hfe,mu,En_F,*SE_L):
	"""	Matsubara sums to be used in n and mu calculation 
	En_F contains negative frequencies only
	SE_L is an optional list of dynamic self-energies """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps] = params_F 
	## define lambdas
	SFband  = lambda x: SFunctionBand(GammaR,GammaL,Delta,x+p.offset_x)
	SFgap   = lambda x: SFunctionGap(GammaR,GammaL,Delta,x+p.offset_x)
	DFband  = lambda x: DeltaFunctionBand(GammaR,GammaL,Delta,Phi,x+p.offset_x)
	DFgap   = lambda x: DeltaFunctionGap(GammaR,GammaL,Delta,Phi,x+p.offset_x)
	Detband = lambda x: DetBand(params_F,hfe,mu,x)
	Detgap  = lambda x: DetGap(params_F,hfe,mu,x)
	## find special points
	dE = sp.around(En_F[1]-En_F[0],8)
	dE_dec = int(-sp.log10(dE))
	EdgePos = sp.nonzero(En_F == sp.around(-Delta,dE_dec))[0][0]
	## fill arrays using lambdas
	Det_F = sp.zeros_like(En_F,dtype = complex)
	SF_F  = sp.zeros_like(En_F,dtype = complex)
	DF_F  = sp.zeros_like(En_F,dtype = complex)
	SF_F[:EdgePos]  = SFband(En_F[:EdgePos])
	SF_F[EdgePos:]  = SFgap(En_F[EdgePos:])
	DF_F[:EdgePos]  = DFband(En_F[:EdgePos])
	DF_F[EdgePos:]  = DFgap(En_F[EdgePos:])
	Det_F[:EdgePos] = Detband(En_F[:EdgePos])
	Det_F[EdgePos:] = Detgap(En_F[EdgePos:])
	## dynamic self-energies
	if len(SE_L) == 4:
		SEn_F     = SE_L[0]
		SEnStar_F = SE_L[1]
		SEa_F     = SE_L[2]
		SEaStar_F = SE_L[3]
		[GFn1_F,GFa1_F,Det1_F] = GreensFunction(params_F,hfe,mu,SEn_F[:EdgePos],SEa_F[:EdgePos],\
			SEnStar_F[:EdgePos],SEaStar_F[:EdgePos],En_F[:EdgePos],'band')
		[GFn2_F,GFa2_F,Det2_F] = GreensFunction(params_F,hfe,mu,SEn_F[EdgePos:],SEa_F[EdgePos:],\
			SEnStar_F[EdgePos:],SEaStar_F[EdgePos:],En_F[EdgePos:],'gap')
		Det_F[:EdgePos] = Det1_F
		Det_F[EdgePos:] = Det2_F
	## Hartree-Fock only
	else:
		SEnStar_F = SEa_F = sp.zeros_like(En_F)
	#for i in range(len(En_F)):
	#	print En_F[i],'\t',sp.real(SF_F[i]),'\t',sp.imag(SF_F[i]),'\t',sp.real(DF_F[i]),'\t',sp.imag(DF_F[i])
	#exit()
	## Matsubara sum of 1/D(iw):
	Int1_F = sp.imag(1.0/Det_F)
	## Matsubara sum of iw(1+s(iw))/D(iw)
	Int2_F = sp.imag((En_F*(1.0+SF_F)+1.0j*GammaNbar-SEnStar_F)/Det_F)
	## Matsubara sum of sp.conj(Delta_Phi(iw))/D(iw)	
	Int3_F = sp.imag((DF_F-SEa_F)/Det_F)
	## Tails
	Tail1 =  Int1_F[0]*En_F[0]/2.0	# behaves as -1/x^3 CHECK!!!!
	Tail2 = -Int2_F[0]*En_F[0]	    # behaves as  1/x^2 CHECK!!!!
	Tail3 = -Int3_F[0]*En_F[0]/2.0	# behaves as  1/x^3 CHECK!!!!
	M1 = -(simps(Int1_F,En_F)+Tail1)/sp.pi
	M2 = -(simps(Int2_F,En_F)+Tail2)/sp.pi
	M3 = -(simps(Int3_F,En_F)+Tail3)/sp.pi
	return [sp.real_if_close(M1),sp.real_if_close(M2),sp.real_if_close(M3)]


def ElectronDensity(params_F,n,mu,En_F,SEn_F,SEa_F):
	"""	calculating n from Matsbara sums MSums """
	N = len(En_F)
	U = params_F[0]
	ed = params_F[6] - U/2.0
	hfe = ed+U*n
	SEnStar_F = -sp.flipud(sp.conj(SEn_F))	# hole self-energy
	SEaStar_F =  sp.flipud(sp.conj(SEa_F))
	MSums_F = MSums(params_F,hfe,mu,En_F[:N/2+1],SEn_F[:N/2+1],SEnStar_F[:N/2+1],\
	SEa_F[:N/2+1],SEaStar_F[:N/2+1])
	n = sp.real_if_close((MSums_F[1]+ed*MSums_F[0])/(1.0-U*MSums_F[0]))
	if sp.imag(n) > 0.0: print '# Warning: non-zero imaginary part \
	of n: {0: .5f}, discarding imag. part.'.format(float(sp.imag(n)))
	return sp.real(n)


def CooperPairDensity(params_F,n,mu,En_F,SEn_F,SEa_F):
	""" calculating mu from Matsbara sums MSums """
	N = len(En_F)
	U = params_F[0]
	ed = params_F[6] - U/2.0
	hfe = ed+U*n
	SEnStar_F = -sp.flipud(sp.conj(SEn_F))	# hole self-energies
	SEaStar_F =  sp.flipud(sp.conj(SEa_F))
	MSums_F = MSums(params_F,hfe,mu,En_F[:N/2+1],SEn_F[:N/2+1],SEnStar_F[:N/2+1],\
	SEa_F[:N/2+1],SEaStar_F[:N/2+1])
	mu =  sp.real_if_close(-MSums_F[2]/(1.0-U*MSums_F[0]))
	if sp.imag(mu) > 0.0: print '# Warning: non-zero imaginary part \
	of mu: {0: .5f}, discarding imag. part.'.format(float(sp.imag(mu)))
	return sp.real(mu)


#####################################################################
# functions to fill arrays ##########################################

def FillGreenHF(params_F,hfe,mu,En_F):
	""" filling the En_F array with HF Green's functions """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps] = params_F 
	## define GFs as lambdas
	GFn_band = lambda x: GFnBand(params_F,hfe,mu,x)
	GFa_band = lambda x: GFaBand(params_F,hfe,mu,x)
	GFn_gap  = lambda x: GFnGap(params_F,hfe,mu,x)
	GFa_gap  = lambda x: GFaGap(params_F,hfe,mu,x)
	## find the frequency mesh
	dE = sp.around(En_F[1]-En_F[0],8)
	dE_dec = int(-sp.log10(dE))
	zero_F = sp.array([0.0])	# zero as an array
	## find special points
	EdgePos1 = sp.nonzero(En_F == sp.around(-Delta,dE_dec))[0][0]
	EdgePos2 = sp.nonzero(En_F == sp.around( Delta,dE_dec))[0][0]
	## fill the arrays
	GFn_F = sp.concatenate((GFn_band(En_F[:EdgePos1]),zero_F,GFn_gap(En_F[EdgePos1+1:EdgePos2])\
	,zero_F,GFn_band(En_F[EdgePos2+1:])))
	GFa_F = sp.concatenate((GFa_band(En_F[:EdgePos1]),zero_F,GFa_gap(En_F[EdgePos1+1:EdgePos2])\
	,zero_F,GFa_band(En_F[EdgePos2+1:])))
	return [GFn_F,GFa_F,EdgePos1,EdgePos2]


def FillGreensFunction(params_F,n,mu,SEn_F,SEa_F,En_F):
	""" calculating the interacting Green's function using the Dyson equation """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps] = params_F 
	hfe = eps+U*(n-0.5)
	dE = sp.around(En_F[1]-En_F[0],8)
	dE_dec = int(-sp.log10(dE))
	EdgePos1 = sp.nonzero(En_F == sp.around(-Delta,dE_dec))[0][0]
	EdgePos2 = sp.nonzero(En_F == sp.around( Delta,dE_dec))[0][0]
	## hole / anti-cooperon self-energies
	SEnStar_F = -sp.flipud(sp.conj(SEn_F))
	SEaStar_F =  sp.flipud(sp.conj(SEa_F))
	## fill the arrays
	[GFn1_F,GFa1_F,Det1_F] = GreensFunction(params_F,hfe,mu,SEn_F[:EdgePos1],SEa_F[:EdgePos1],\
		SEnStar_F[:EdgePos1],SEaStar_F[:EdgePos1],En_F[:EdgePos1],'band')
	[GFn2_F,GFa2_F,Det2_F] = GreensFunction(params_F,hfe,mu,SEn_F[EdgePos1+1:EdgePos2],SEa_F[EdgePos1+1:EdgePos2],\
		SEnStar_F[EdgePos1+1:EdgePos2],SEaStar_F[EdgePos1+1:EdgePos2],En_F[EdgePos1+1:EdgePos2],'gap')
	[GFn3_F,GFa3_F,Det3_F] = GreensFunction(params_F,hfe,mu,SEn_F[EdgePos2+1:],SEa_F[EdgePos2+1:],\
		SEnStar_F[EdgePos2+1:],SEaStar_F[EdgePos2+1:],En_F[EdgePos2+1:],'band')
	## concatenate the band and gap parts together
	Det_F = sp.concatenate([Det1_F,sp.zeros(1),Det2_F,sp.zeros(1),Det3_F])
	GFn_F = sp.concatenate([GFn1_F,sp.zeros(1),GFn2_F,sp.zeros(1),GFn3_F])
	GFa_F = sp.concatenate([GFa1_F,sp.zeros(1),GFa2_F,sp.zeros(1),GFa3_F])
	return [GFn_F,GFa_F,Det_F]


#####################################################################
# Auxiliary functions ###############################################

def PrintDet(params_F,hfe,mu,En_F):
	""" print the determinant for testing purposes """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps] = params_F 
	dE = sp.around(En_F[1]-En_F[0],8)
	dE_dec = int(-sp.log10(dE))
	EdgePos1 = sp.nonzero(En_F == sp.around(-Delta,dE_dec))[0][0]
	EdgePos2 = sp.nonzero(En_F == sp.around(Delta,dE_dec))[0][0]
	Det_F = sp.zeros_like(En_F,dtype = complex)
	Detband = lambda x: DetBand(params_F,hfe,mu,x)
	Detgap  = lambda x: DetGap(params_F,hfe,mu,x)
	Det_F[:EdgePos1] = Detband(En_F[:EdgePos1])
	Det_F[EdgePos1:EdgePos2] = Detgap(En_F[EdgePos1:EdgePos2])
	Det_F[EdgePos2:] = Detband(En_F[EdgePos2:])
	for i in range(len(En_F)):
		print En_F[i],'\t',sp.real(Det_F[i]),'\t',sp.imag(Det_F[i])

## ssnlib.py end ##

