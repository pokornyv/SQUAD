# SQUAD - superconducting quantum dot
# functions library for case with two sc and one normal lead
# uses scipy, optimized on Python 2.7.5
# Vladislav Pokorny, 2016; pokornyv@fzu.cz

from sys import version_info,argv
from time import ctime
import scipy as sp
from scipy.integrate import simps
from scipy.fftpack import fft,ifft
from scipy.interpolate import InterpolatedUnivariateSpline
from squadlib1 import SFunctionBand,SFunctionGap,DeltaFunctionBand,DeltaFunctionGap
from squadlib1 import FermiDirac,BoseEinstein,FindEdges,FindInEnergies,Brillouin
from squadlib2 import KramersKronigFFT
import params as p

#####################################################################
# Green's function determinants on real axis ########################

def DetBand(params_F,hfe,hm,mu,s,x):
	"""	determinant of the HF Green's function in 
	the band region (-inf:-Delta),(Delta,inf) """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	SF = lambda x: SFunctionBand(GammaR,GammaL,Delta,x+p.P['offset_x'])
	DF = lambda x: DeltaFunctionBand(GammaR,GammaL,Delta,Phi,x+p.P['offset_x'])
	return (x*(1.0+SF(x))+s*h+1.0j*GammaNbar)**2-hfe**2-(DF(x)-U*mu)**2


def DetGap(params_F,hfe,hm,mu,s,x):
	"""	determinant of the HF Green's function in 
	the gap region (-Delta:Delta) """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	SF = lambda x: SFunctionGap(GammaR,GammaL,Delta,x+p.P['offset_x'])
	DF = lambda x: DeltaFunctionGap(GammaR,GammaL,Delta,Phi,x+p.P['offset_x'])
	return (x*(1.0+SF(x))+s*hm+1.0j*GammaNbar)**2-hfe**2-(DF(x)-U*mu)**2 


def GFnBand(params_F,hfe,hm,mu,s,x):
	""" normal GF in band region """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	SF = lambda x: SFunctionBand(GammaR,GammaL,Delta,x+p.P['offset_x'])
	DF = lambda x: DeltaFunctionBand(GammaR,GammaL,Delta,Phi,x+p.P['offset_x'])
	Det = lambda x: DetBand(params_F,hfe,hm,mu,s,x)
	return (x*(1.0+SF(x))+s*hm+1.0j*GammaNbar+hfe)/Det(x)


def GFnGap(params_F,hfe,hm,mu,s,x):
	""" normal GF in gap region """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	SF = lambda x: SFunctionGap(GammaR,GammaL,Delta,x+p.P['offset_x'])
	DF = lambda x: DeltaFunctionGap(GammaR,GammaL,Delta,Phi,x+p.P['offset_x'])
	Det = lambda x: DetGap(params_F,hfe,hm,mu,s,x)
	return (x*(1.0+SF(x))+s*hm+1.0j*GammaNbar+hfe)/Det(x)


def GFaBand(params_F,hfe,hm,mu,s,x):
	""" anomalous GF in band region """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	SF = lambda x: SFunctionBand(GammaR,GammaL,Delta,x+p.P['offset_x'])
	DF = lambda x: DeltaFunctionBand(GammaR,GammaL,Delta,Phi,x+p.P['offset_x'])
	Det = lambda x: DetBand(params_F,hfe,hm,mu,s,x)
	return -(DF(x)-U*mu)/Det(x)


def GFaGap(params_F,hfe,hm,mu,s,x):
	""" anomalous GF in gap region """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	SF = lambda x: SFunctionGap(GammaR,GammaL,Delta,x+p.P['offset_x'])
	DF = lambda x: DeltaFunctionGap(GammaR,GammaL,Delta,Phi,x+p.P['offset_x'])
	Det = lambda x: DetGap(params_F,hfe,hm,mu,s,x)
	return -(DF(x)-U*mu)/Det(x)


def GreensFunction(params_F,hfe,hm,mu,s,SEn_F,SEa_F,SEnStar_F,SEaStar_F,En_F,A):
	"""	general form the Nambu Green's function and its determinant
	A='band'/'gap' calculates determinant for band/gap """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	if A == 'band':
		S_F = SFunctionBand(GammaR,GammaL,Delta,En_F+p.P['offset_x'])
		D_F = DeltaFunctionBand(GammaR,GammaL,Delta,Phi,En_F+p.P['offset_x'])
	else:
		S_F = SFunctionGap(GammaR,GammaL,Delta,En_F+p.P['offset_x'])
		D_F = DeltaFunctionGap(GammaR,GammaL,Delta,Phi,En_F+p.P['offset_x'])
	Detn_F = (En_F*(1.0+S_F)+s*hm+1.0j*GammaNbar-hfe-SEn_F)*(En_F*(1.0+S_F)+s*hm+1.0j*GammaNbar+hfe-SEnStar_F)
	Deta_F = (D_F-U*mu-SEa_F)*(D_F-U*mu-SEaStar_F)
	Det_F = Detn_F-Deta_F
	GFn_F =  (En_F*(1.0+S_F)+s*hm+1.0j*GammaNbar+hfe-SEnStar_F)/Det_F
	GFa_F = -(D_F-U*mu-SEa_F)/Det_F
	return [GFn_F,GFa_F,Det_F]


#####################################################################
# Hartree-Fock solver ###############################################

def SolveHFssn(params_F):
	""" HF equations solver for case with two SC and one normal lead """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	s  = 1 
	m = 0.0 			# for compatibility with magnetic solution
	hm = h+U*m/2.0
	ed = eps-U/2.0		# localized energy level shifted to symmetry point
	ErrMsg = 0			# error message indicator
	## fill the energy array
	dE   = 1e-4							# band energy sampling
	Emin = -100.0						# lower cutoff for band energy
	En_F = sp.around(sp.arange(Emin,0.0,dE),int(-sp.log10(dE)))
	## initial condition: normal state w/o superconductivity
	from scipy.optimize import fixed_point
	Gamma = (GammaR + GammaL + GammaNbar)
	n_density = lambda x: 0.5 - sp.arctan((ed+U*x)/Gamma)/sp.pi
	n = fixed_point(n_density,0.5)
	mu = 0.0 if Delta == 0 else 0.2	# good guess in zero-phase, pi-phase needs small and negative
	hfe = ed+U*n
	if Delta != 0:
		n_old = 1e5
		mu_old = 1e5
		i = 0		# iterations counter
		imax = p.P['HF_max_iter']	# maximum number of iterations
		while any([sp.fabs(mu-mu_old) > p.P['ConvHF'], sp.fabs(n-n_old) > p.P['ConvHF']]):
			n_old = n
			mu_old = mu
			hfe = ed+U*n
			[D1,D2,D3] = MSums(params_F,hfe,hm,mu,s,En_F)
			mu = -D3/(1.0-U*D1)
			[D1,D2,D3] = MSums(params_F,hfe,hm,mu,s,En_F)
			if eps == 0: n = 0.5
			else: n = (D2+ed*D1)/(1.0-U*D1)
			hfe = ed+U*n
			if i > imax: 
				print('# Warning: SolveHF: No convergence after '+str(i)+' iterations, exit.')
				n = mu = -1.0
				ErrMsg = 1	
				break
			if p.P['WriteIO'] and i+1%100 == 0: print('# {0: 3d}\t n = {1: .6f} +{2: .6f}i,  mu = {3: .6f} +{4: .6f}i'\
			.format(i+1,float(sp.real(n)),float(sp.imag(n)),float(sp.real(mu)),float(sp.imag(mu))))
			#if WriteIO: print('# {0: 3d}\t n = {1: .6f} +{2: .6f}i,  mu = {3: .6f} +{4: .6f}i'\
			#.format(i+1,float(sp.real(n)),float(sp.imag(n)),float(sp.real(mu)),float(sp.imag(mu))))
			i += 1
		if p.P['WriteIO']: print('# {0: 3d} iterations,\t n = {1: .6f} +{2: .6f}i,  mu = {3: .6f} +{4: .6f}i'\
		.format(i,float(sp.real(n)),float(sp.imag(n)),float(sp.real(mu)),float(sp.imag(mu))))
	return sp.array([n,mu,ErrMsg])


def SolveHFssnMag(params_F):
	""" HF equations solver for case with two SC and one normal lead 
	including magnetic field. Should replace SolveHFssn after testing"""
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	ed = eps-U/2.0		# localized energy level shifted to symmetry point
	ErrMsg = 0			# error message indicator
	## fill the energy array
	dE   = 1e-4							# band energy sampling
	Emin = -100.0						# lower cutoff for band energy
	En_F = sp.around(sp.arange(Emin,0.0,dE),int(-sp.log10(dE)))
	## initial condition: normal state w/o superconductivity
	from scipy.optimize import fixed_point
	Gamma = (GammaR + GammaL + GammaNbar)
	n = 0.5
	m = Brillouin(0.5,h)
	[n_old,m_old] = [1e5,1e5]
	while any([sp.fabs(m_old-m) > p.P['ConvHF'],sp.fabs(n_old-n) > p.P['ConvHF']]):
		n_old = n
		m_old = m
		nf = lambda x: 0.5*(1.0 - (sp.arctan((ed+h+U*(x+m))/Gamma) + sp.arctan((ed-h+U*(x-m))/Gamma))/sp.pi)
		n = fixed_point(nf,0.5)
		if h!=0.0:
			mf = lambda x: (sp.arctan((ed+h+U*(n+x))/Gamma) - sp.arctan((ed-h+U*(n-x))/Gamma))/sp.pi
			m = fixed_point(mf,0.99*sp.sign(h))
	mu = 0.0 if Delta == 0 else 0.2	# good guess in zero-phase, pi-phase needs small and negative
	print 'n0,m0,mu0: ',m,n,mu
	## calculate the superconducting case #################
	if Delta != 0:
		[n_old,m_old,mu_old] = [1e5,1e5,1e5]
		i = 0		# iterations counter
		imax = p.P['HF_max_iter']	# maximum number of iterations
		while any([sp.fabs(mu-mu_old) > p.P['ConvHF'], \
			sp.fabs(n-n_old) > p.P['ConvHF'], \
			sp.fabs(m-m_old) > p.P['ConvHF']]):
			[n_old,m_old,mu_old] = [n,m,mu]
			hfe = ed+U*n
			hm  = h+U*m/2.0
			[D1up,D2up,D3up] = MSums(params_F,hfe,hm,mu,1,En_F)
			[D1dn,D2dn,D3dn] = [D1up,D2up,D3up] if h == 0 else MSums(params_F,hfe,hm,mu,-1,En_F)
			mu_up = -D3up/(1.0-U*D1up)
			mu_dn = -D3dn/(1.0-U*D1dn)
			print 'mus:',mu_up,mu_dn
			# mu_up and mu_dn must be the same, next line is just to smooth out possible numerical errors
			mu = 0.5*(mu_up+mu_dn)
			[D1up,D2up,D3up] = MSums(params_F,hfe,hm,mu,1,En_F)
			[D1dn,D2dn,D3dn] = [D1up,D2up,D3up] if h == 0 else MSums(params_F,hfe,hm,mu,-1,En_F)
			n = 0.5 if eps == 0.0 else 0.5*(D2up+D2dn+D1up*(ed+hm)+D1dn*(ed-hm))/(1.0-U/2.0*(D1up+D1dn))
			m = 0.0	if h   == 0.0 else D2up-D2dn+D1up*(hfe+h)-D1dn*(hfe-h)/(1.0-U/2.0*(D1up+D1dn))
			if i > imax: 
				print('# Warning: SolveHF: No convergence after '+str(i)+' iterations, exit.')
				n = mu = m = -1.0
				ErrMsg = 1	
				break
			print i,n,m,mu
			if p.P['WriteIO'] and i+1%10 == 0: 
				print('# {0: 3d}\t n = {1: .6f} +{2: .6f}i,  m = {3: .6f} +{4: .6f}i,  mu = {5: .6f} +{6: .6f}i'\
				.format(i+1,float(sp.real(n)),float(sp.imag(n)),float(sp.real(m)),float(sp.imag(m)),\
				float(sp.real(mu)),float(sp.imag(mu))))
			i += 1
		if p.P['WriteIO']: \
		print('# {0: 3d} iterations,\tn = {1: .6f} +{2: .6f}i,\tm = {3: .6f} +{4: .6f}i,\tmu = {5: .6f} +{6: .6f}i'\
		.format(i,float(sp.real(n)),float(sp.imag(n)),\
		float(sp.real(m)),float(sp.imag(m)),float(sp.real(mu)),float(sp.imag(mu))))
	return sp.array([n,m,mu,ErrMsg])


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

def MSums(params_F,hfe,hm,mu,s,En_F,*SE_L):
	"""	Matsubara sums to be used in n and mu calculation 
	En_F contains negative frequencies only
	SE_L is an optional list of dynamic self-energies """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	## define lambdas
	SFband  = lambda x: SFunctionBand(GammaR,GammaL,Delta,x+p.P['offset_x'])
	SFgap   = lambda x: SFunctionGap(GammaR,GammaL,Delta,x+p.P['offset_x'])
	DFband  = lambda x: DeltaFunctionBand(GammaR,GammaL,Delta,Phi,x+p.P['offset_x'])
	DFgap   = lambda x: DeltaFunctionGap(GammaR,GammaL,Delta,Phi,x+p.P['offset_x'])
	Detband = lambda x: DetBand(params_F,hfe,hm,mu,s,x)
	Detgap  = lambda x: DetGap(params_F,hfe,hm,mu,s,x)
	## find special points
	dE = sp.around(En_F[1]-En_F[0],8)
	dE_dec = int(-sp.log10(dE))
	if Delta < sp.fabs(En_F[0]): EdgePos = sp.nonzero(En_F == sp.around(-Delta,dE_dec))[0][0]
	else: EdgePos = 0
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
	if len(SE_L) == 4: 	## dynamic self-energies
		SEn_F     = SE_L[0]
		SEnStar_F = SE_L[1]
		SEa_F     = SE_L[2]
		SEaStar_F = SE_L[3]
		[GFn1_F,GFa1_F,Det1_F] = GreensFunction(params_F,hfe,hm,mu,s,SEn_F[:EdgePos],SEa_F[:EdgePos],\
			SEnStar_F[:EdgePos],SEaStar_F[:EdgePos],En_F[:EdgePos],'band')
		[GFn2_F,GFa2_F,Det2_F] = GreensFunction(params_F,hfe,hm,mu,s,SEn_F[EdgePos:],SEa_F[EdgePos:],\
			SEnStar_F[EdgePos:],SEaStar_F[EdgePos:],En_F[EdgePos:],'gap')
		Det_F[:EdgePos] = Det1_F
		Det_F[EdgePos:] = Det2_F
	else: 	## Hartree-Fock only
		SEnStar_F = SEa_F = sp.zeros_like(En_F)
	## Matsubara sum of 1/D(iw):
	Int1_F = sp.imag(1.0/Det_F)
	## Matsubara sum of iw(1+s(iw))/D(iw)
	Int2_F = sp.imag((En_F*(1.0+SF_F)+1.0j*GammaNbar-SEnStar_F)/Det_F)
	## Matsubara sum of Delta_Phi(iw)/D(iw)	
	Int3_F = sp.imag((DF_F-SEa_F)/Det_F)
	## Tails
	Tail1 =  Int1_F[0]*En_F[0]/2.0	# behaves as -1/x^3 CHECK!!!!
	Tail2 = -Int2_F[0]*En_F[0]	    # behaves as  1/x^2 CHECK!!!!
	Tail3 = -Int3_F[0]*En_F[0]/2.0	# behaves as  1/x^3 CHECK!!!!
	M1 = -(simps(Int1_F,En_F)+Tail1)/sp.pi
	M2 = -(simps(Int2_F,En_F)+Tail2)/sp.pi
	M3 = -(simps(Int3_F,En_F)+Tail3)/sp.pi
	return [sp.real_if_close(M1),sp.real_if_close(M2),sp.real_if_close(M3)]


def ElectronDensity(params_F,n,m,mu,s,En_F,SEn_F,SEa_F):
	"""	calculating n from Matsbara sums MSums """
	N = len(En_F)
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	ed  = eps-U/2.0
	hfe = ed+U*n
	hm  = h+U*m/2.0
	SEnStar_F = -sp.flipud(sp.conj(SEn_F))	# hole self-energy
	SEaStar_F =  sp.flipud(sp.conj(SEa_F))
	MSums_F = MSums(params_F,hfe,hm,mu,s,En_F[:N/2+1],SEn_F[:N/2+1],SEnStar_F[:N/2+1],\
	SEa_F[:N/2+1],SEaStar_F[:N/2+1])
	n = sp.real_if_close((MSums_F[1]+ed*MSums_F[0])/(1.0-U*MSums_F[0]))
	if sp.imag(n) > 0.0: print('# Warning: non-zero imaginary part \
	of n: {0: .5f}, discarding imag. part.'.format(float(sp.imag(n))))
	return sp.real(n)


def CooperPairDensity(params_F,n,m,mu,s,En_F,SEn_F,SEa_F):
	""" calculating mu from Matsbara sums MSums """
	N = len(En_F)
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	ed  = eps-U/2.0
	hfe = ed+U*n
	hm  = h+U*m/2.0
	SEnStar_F = -sp.flipud(sp.conj(SEn_F))	# hole self-energies
	SEaStar_F =  sp.flipud(sp.conj(SEa_F))
	MSums_F = MSums(params_F,hfe,hm,mu,s,En_F[:N/2+1],SEn_F[:N/2+1],SEnStar_F[:N/2+1],\
	SEa_F[:N/2+1],SEaStar_F[:N/2+1])
	mu =  sp.real_if_close(-MSums_F[2]/(1.0-U*MSums_F[0]))
	if sp.imag(mu) > 0.0: print('# Warning: non-zero imaginary part \
	of mu: {0: .5f}, discarding imag. part.'.format(float(sp.imag(mu))))
	return sp.real(mu)


#####################################################################
# Josephson current and conductance #################################

def JosephsonCurrent(GFa_F,En_F,params_F,direction):
	""" Josephson current flowing from 'direction' lead to the other one
	direction = 'R' (right) or 'L' (left) 
	returns separately the band and the gap contributions """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	Gamma = GammaR if direction == 'R' else -GammaL
	MultFac = -4.0*Delta*Gamma*sp.sin(Phi/2.0)
	## define lambdas
	SFband  = lambda x: SFunctionBand(GammaR,GammaL,Delta,x+p.P['offset_x'])
	SFgap   = lambda x: SFunctionGap(GammaR,GammaL,Delta,x+p.P['offset_x'])
	# find band edge
	[EdgePos1,EdgePos2] = FindEdges(En_F,Delta)
	# integrate separately band and gap
	if Delta < sp.fabs(En_F[0]): # finite gap
		Int1_F  = sp.real(GFa_F[:EdgePos1])*sp.imag(SFband(En_F[:EdgePos1]))
		I1 = sp.real_if_close(simps(Int1_F,En_F[:EdgePos1]))/sp.pi
	else: I1 = 0.0
	Int2_F  = sp.imag(GFa_F[EdgePos1:len(En_F)/2+1])*SFgap(En_F[EdgePos1:len(En_F)/2+1])
	I2 = sp.real_if_close(simps(Int2_F,En_F[EdgePos1:len(En_F)/2+1]))/sp.pi
	return [MultFac*I1,MultFac*I2]
	

def AndreevConductance(w,GFa_F,En_F,params_F):
	""" Andreev conductance between SC lead and normal lead """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F
	GammaN = 2.0*GammaNbar
	Pos = FindInEnergies(En_F,w)
	GFa = GFa_F[Pos]
	return sp.real_if_close(4.0*GammaN**2*GFa*sp.conj(GFa))


def WriteAndreevConductance(GFa_F,GFn_F,En_F,params_F):
	""" differential conductance between SC lead and normal lead 
	calculated for all subgap values second column with normal spectral
	function is added so we can compare the peak widths """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F
	GammaN = 2.0*GammaNbar # GammaN, not GammaN/2 enters the calculation
	D = p.P['EmaxFiles'] if Delta > p.P['EmaxFiles'] else Delta # taking care of infinite Delta
	# find band edge and select gap region only
	[EdgePos1,EdgePos2] = FindEdges(En_F,D)
	Gap_F = En_F[EdgePos1:EdgePos2+1]
	GapGFa_F = GFa_F[EdgePos1:EdgePos2+1]
	GapDOS_F = -sp.imag(GFn_F[EdgePos1:EdgePos2+1])/sp.pi
	# calculate the conductance
	AC_F = sp.real_if_close(4.0*GammaN**2*GapGFa_F*sp.conj(GapGFa_F))
	f = open('AC_P'+str(Phi/sp.pi)+'U'+str(U)+'.dat','w')
	ver = str(version_info[0])+'.'+str(version_info[1])+'.'+str(version_info[2])
	text_header = '# generated by '+str(argv[0])+', python version: '+str(ver)\
	+', SciPy version: '+str(sp.version.version)+'\n# '+ctime()+'\n'
	f.write(text_header)
	# write conductance to file
	for i in range(len(Gap_F)):
		if i % p.P['EstepFiles'] == 0 or sp.fabs(Gap_F[i]) < 1e-6:
			f.write('{0: .4f}\t{1: .8f}\t{2: .8f}\n'\
			.format(float(Gap_F[i]),float(AC_F[i]),float(GapDOS_F[i])))
	f.close()
	return AC_F


#####################################################################
# functions to fill arrays ##########################################

def FillGreenHF(params_F,hfe,hm,mu,s,En_F):
	""" filling the En_F array with HF Green's functions """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	## define GFs as lambdas
	GFn_band = lambda x: GFnBand(params_F,hfe,hm,mu,s,x)
	GFa_band = lambda x: GFaBand(params_F,hfe,hm,mu,s,x)
	GFn_gap  = lambda x: GFnGap(params_F,hfe,hm,mu,s,x)
	GFa_gap  = lambda x: GFaGap(params_F,hfe,hm,mu,s,x)
	## find the frequency mesh
	zero_F = sp.array([0.0])	# zero as an array
	## find special points
	[EdgePos1,EdgePos2] = FindEdges(En_F,Delta)
	## fill the arrays
	GFn_F = sp.concatenate((GFn_band(En_F[:EdgePos1]),zero_F,GFn_gap(En_F[EdgePos1+1:EdgePos2])\
	,zero_F,GFn_band(En_F[EdgePos2+1:])))
	GFa_F = sp.concatenate((GFa_band(En_F[:EdgePos1]),zero_F,GFa_gap(En_F[EdgePos1+1:EdgePos2])\
	,zero_F,GFa_band(En_F[EdgePos2+1:])))
	return [GFn_F,GFa_F]


def FillGreensFunction(params_F,n,m,mu,s,SEn_F,SEa_F,En_F):
	""" calculating the interacting Green's function using the Dyson equation """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	hfe = eps+U*(n-0.5)
	hm = h+U*m/2.0
	[EdgePos1,EdgePos2] = FindEdges(En_F,Delta)
	## hole / anti-cooperon self-energies
	SEnStar_F = -sp.flipud(sp.conj(SEn_F))
	SEaStar_F =  sp.flipud(sp.conj(SEa_F))
	## fill the arrays
	[GFn1_F,GFa1_F,Det1_F] = GreensFunction(params_F,hfe,hm,mu,s,SEn_F[:EdgePos1],SEa_F[:EdgePos1],\
		SEnStar_F[:EdgePos1],SEaStar_F[:EdgePos1],En_F[:EdgePos1],'band')
	[GFn2_F,GFa2_F,Det2_F] = GreensFunction(params_F,hfe,hm,mu,s,SEn_F[EdgePos1+1:EdgePos2],SEa_F[EdgePos1+1:EdgePos2],\
		SEnStar_F[EdgePos1+1:EdgePos2],SEaStar_F[EdgePos1+1:EdgePos2],En_F[EdgePos1+1:EdgePos2],'gap')
	[GFn3_F,GFa3_F,Det3_F] = GreensFunction(params_F,hfe,hm,mu,s,SEn_F[EdgePos2+1:],SEa_F[EdgePos2+1:],\
		SEnStar_F[EdgePos2+1:],SEaStar_F[EdgePos2+1:],En_F[EdgePos2+1:],'band')
	## concatenate the band and gap parts together
	Det_F = sp.concatenate([Det1_F,sp.zeros(1),Det2_F,sp.zeros(1),Det3_F])
	GFn_F = sp.concatenate([GFn1_F,sp.zeros(1),GFn2_F,sp.zeros(1),GFn3_F])
	GFa_F = sp.concatenate([GFa1_F,sp.zeros(1),GFa2_F,sp.zeros(1),GFa3_F])
	return [GFn_F,GFa_F,Det_F]


#####################################################################
# auxiliary functions ###############################################

def PrintDet(params_F,hfe,hm,mu,En_F):
	""" print the determinant for testing purposes """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	[EdgePos1,EdgePos2] = FindEdges(En_F,Delta)
	Det_F = sp.zeros_like(En_F,dtype = complex)
	Detband = lambda x: DetBand(params_F,hfe,hm,mu,x)
	Detgap  = lambda x: DetGap(params_F,hfe,hm,mu,x)
	Det_F[:EdgePos1] = Detband(En_F[:EdgePos1])
	Det_F[EdgePos1:EdgePos2] = Detgap(En_F[EdgePos1:EdgePos2])
	Det_F[EdgePos2:] = Detband(En_F[EdgePos2:])
	for i in range(len(En_F)):
		print('{0: .6f}\t{1: .6f}\t{2: .6f}'\
		.format(float(En_F[i]),float(sp.real(Det_F[i])),float(sp.imag(Det_F[i]))))


def PrintHybs(params_F,En_F):
	""" prints the hybridization functions S and Delta """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	[EdgePos1,EdgePos2] = FindEdges(En_F,Delta)
	S1_F = lambda x: SFunctionBand(GammaR,GammaL,Delta,x+p.P['offset_x'])
	D1_F = lambda x: DeltaFunctionBand(GammaR,GammaL,Delta,Phi,x+p.P['offset_x'])
	S2_F = lambda x: SFunctionGap(GammaR,GammaL,Delta,x+p.P['offset_x'])
	D2_F = lambda x: DeltaFunctionGap(GammaR,GammaL,Delta,Phi,x+p.P['offset_x'])
	for w in En_F[:EdgePos1]:
		print('{0: .4f}\t{1: .8f}\t{2: .8f}\t{3: .8f}\t{4: .8f}')\
		.format(float(w),float(sp.real(S1_F(w))),float(sp.imag(S1_F(w)))\
		,float(sp.real(D1_F(w))),float(sp.imag(D1_F(w))))
	for w in En_F[EdgePos1:EdgePos2+1]:
		print('{0: .4f}\t{1: .8f}\t{2: .8f}\t{3: .8f}\t{4: .8f}')\
		.format(float(w),float(sp.real(S2_F(w))),float(sp.imag(S2_F(w)))\
		,float(sp.real(D2_F(w))),float(sp.imag(D2_F(w))))
	for w in En_F[EdgePos2+1:]:
		print('{0: .4f}\t{1: .8f}\t{2: .8f}\t{3: .8f}\t{4: .8f}')\
		.format(float(w),float(sp.real(S1_F(w))),float(sp.imag(S1_F(w)))\
		,float(sp.real(D1_F(w))),float(sp.imag(D1_F(w))))


def ZeroEnergyDiff(GFn_F,En_F,params_F):
	""" determines if the value at Fermi energy is a minimum or maximum """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F
	[DPos1,DPos2] = FindEdges(En_F,Delta/4.0)
	F = InterpolatedUnivariateSpline(En_F[DPos1:DPos2],-sp.imag(GFn_F[DPos1:DPos2])/sp.pi)
	return [F(0),F.derivatives(0)[1],F.derivatives(0)[2]]


def SubgapStatesMaxima(GFn_F,En_F,params_F,WriteIO):
	""" searches for the maxima in the gap region, gives a 
	rough guess of the ABS/Hubbard band positions """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F 
	[EdgePos1,EdgePos2] = FindEdges(En_F,Delta)
	Gap_F = En_F[EdgePos1:EdgePos2]
	GapDoS_F = -sp.imag(GFn_F[EdgePos1:EdgePos2])/sp.pi
	GapDosdiff1_F = sp.zeros_like(GapDoS_F)
	GapDosdiff2_F = sp.zeros_like(GapDoS_F)
	F = InterpolatedUnivariateSpline(Gap_F,GapDoS_F)
	SGapStates_F = sp.empty(0)
	for i in range(len(Gap_F)):
		[GapDosdiff1_F[i],GapDosdiff2_F[i]] = F.derivatives(Gap_F[i])[1:3]
		if all([i>5,i < len(Gap_F)-5,GapDosdiff1_F[i-1]*GapDosdiff1_F[i]<0.0,GapDosdiff2_F[i]<0.0]):
			#print Gap_F[i],GapDoS_F[i],GapDosdiff1_F[i],GapDosdiff2_F[i]
			pos = float((Gap_F[i-1]+Gap_F[i])/2.0)
			if WriteIO: print('# - maximum in spectral function at {0: .4f}'.format(pos))
			SGapStates_F = sp.concatenate([SGapStates_F,sp.array([pos])])
	if WriteIO and len(SGapStates_F) == 0:
		print('# no subgap maxima in spectral function')
	return SGapStates_F

## ssnlib.py end ##

