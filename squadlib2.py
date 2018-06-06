## SQUAD - superconducting quantum dot
## functions library 2
## Vladislav Pokorny; 2012-2018; pokornyv@fzu.cz

from __future__ import print_function
import scipy as sp
from scipy.integrate import trapz,simps
from scipy.fftpack import fft,ifft
from scipy.interpolate import InterpolatedUnivariateSpline,UnivariateSpline
from squadlib1 import SFunctionBand,SFunctionGap,DeltaFunctionBand,DeltaFunctionGap,FindEdges
from sys import version_info,argv,exit
from time import ctime


#####################################################################
# output functions ##################################################


def WriteFile(En_F,Xn_F,Xa_F,params_F,pole_pos,f_type,Emax,NE):
	"""	writes an output file suitable for gnuplot
	range for output is (-Emax:Emax) with step NE x dE
	pole_pos guarantees we don't miss the poles """
	[U,Delta,GammaR,GammaL,GammaNbar,Phi,eps,h] = params_F
	P = sp.around(Phi/sp.pi,3)
	GammaLR = GammaL/GammaR if GammaR != 0.0 else 1.0
	GammaN = 2.0*GammaNbar  ## compatibility with the SSN branch
	[kmin,kmax] = FindEdges(En_F,Emax)
	filename = ""
	## long filenames contain calculation parameters, uncomment if needed (e.g. multiple cases in one folder
	#filename = 'U'+str(U)+'D'+str(Delta)+'GR'+str(GammaR)+'GLR'+str(sp.around(GammaLR,2))+'e'+str(eps)+'P'+str(P)
	f = open(filename+f_type+'.dat','w')
	text_header = '# U='+str(U)+', Delta='+str(Delta)+', eps='+str(eps)+', Phi/pi='+str(P)\
	+'\n# GammaR='+str(GammaR)+', GammaL='+str(GammaL)+', GammaL/GammaR='+str(GammaLR)+', GammaN='+str(GammaN)+'\n'
	f.write(text_header)
	ver = str(version_info[0])+'.'+str(version_info[1])+'.'+str(version_info[2])
	text_header = '# generated by '+str(argv[0])+', python version: '+str(ver)\
	+', SciPy version: '+str(sp.version.version)+'\n# '+ctime()+'\n'
	f.write(text_header)
	[xzeroPos1,xzeroPos2] = FindEdges(En_F,pole_pos)
	for k in range(len(En_F[kmin:kmax])):
		if any([k%NE == 0,k == xzeroPos1-kmin,k == xzeroPos2-kmin]):
			gline = '{0:.5f}\t{1:.8f}\t{2:.8f}\t{3:.8f}\t{4:.8f}'\
			.format(sp.float32(En_F[k+kmin]),sp.float32(sp.real(Xn_F[k+kmin])),sp.float32(sp.imag(Xn_F[k+kmin]))\
			,sp.float32(sp.real(Xa_F[k+kmin])),sp.float32(sp.imag(Xa_F[k+kmin])))+'\n'
			f.write(gline)	
	f.close()


#####################################################################
# convolution procedures using FFT ##################################

def TwoParticleBubbles(GFn_F,GFa_F,Zn_F,wzero):
	""" calcualting two-particle bubbles,
	 uses fft from scipy.fftpack """
	N = (len(Zn_F)-1)/2
	dz = sp.around(Zn_F[1]-Zn_F[0],8)
	dz_dec = int(-sp.log10(dz))
	GFnh_F = -sp.conj(sp.flipud(GFn_F))	
	## zero-padding of the arrays to double the size
	FD_F = 1.0*sp.concatenate([sp.zeros(3*N+3),sp.ones(N+1)])
	ImGFp_F = sp.concatenate([sp.imag(GFn_F[(len(Zn_F)-1)/2:])\
	,sp.zeros(2*N+3),sp.imag(GFn_F[:(len(Zn_F)-1)/2])])
	ImGFh_F = sp.concatenate([sp.imag(GFnh_F[(len(Zn_F)-1)/2:])\
	,sp.zeros(2*N+3),sp.imag(GFnh_F[:(len(Zn_F)-1)/2])])
	ImGFa_F = sp.concatenate([sp.imag(GFa_F[(len(Zn_F)-1)/2:])\
	,sp.zeros(2*N+3),sp.imag(GFa_F[:(len(Zn_F)-1)/2])])
	## perform convolution/cross-correlation via FFT 
	ftImChin1_F = -sp.conj(fft(FD_F*ImGFp_F))*fft(ImGFp_F)*dz
	ftImChin2_F = fft(FD_F*ImGFp_F)*sp.conj(fft(ImGFp_F))*dz
	ImChin_F = sp.real(ifft(ftImChin1_F+ftImChin2_F))/sp.pi
	ImChin_F = sp.concatenate([ImChin_F[3*N+4:],ImChin_F[:N+1]])
	ftImChia1_F = -sp.conj(fft(FD_F*ImGFa_F))*fft(ImGFa_F)*dz
	ftImChia2_F = fft(FD_F*ImGFa_F)*sp.conj(fft(ImGFa_F))*dz
	ImChia_F = sp.real(ifft(ftImChia1_F+ftImChia2_F))/sp.pi
	ImChia_F = sp.concatenate([ImChia_F[3*N+4:],ImChia_F[:N+1]])
	## find ABS positions 2w_0 on the energy axis
	if sp.fabs(wzero)>dz:
		ABSposChi1 = sp.nonzero(Zn_F == -2.0*sp.around(wzero,dz_dec))[0][0]
		ABSposChi2 = sp.nonzero(Zn_F ==  2.0*sp.around(wzero,dz_dec))[0][0]
	else:	## pi-phase Green's functions, putting poles at lowest possible points
		ABSposChi1 = sp.nonzero(Zn_F == -2.0*sp.around(dz,dz_dec))[0][0]
		ABSposChi2 = sp.nonzero(Zn_F ==  2.0*sp.around(dz,dz_dec))[0][0]	
	## add residues to given positions
	ResChin1 = -ImChin_F[ABSposChi1]*dz/sp.pi
	ResChin2 = -ImChin_F[ABSposChi2]*dz/sp.pi
	ResChia1 = -ImChia_F[ABSposChi1]*dz/sp.pi
	ResChia2 = -ImChia_F[ABSposChi2]*dz/sp.pi
	## find real part from imaginary using KK relations
	Chin_F = KramersKronigFFT_ABS(ImChin_F,Zn_F,[ABSposChi1,ABSposChi2],[ResChin1,ResChin2]) + 1.0j*ImChin_F
	Chia_F = KramersKronigFFT_ABS(ImChia_F,Zn_F,[ABSposChi1,ABSposChi2],[ResChia1,ResChia2]) + 1.0j*ImChia_F
	Chin_F[ABSposChi1] = 1.0j*sp.imag(Chin_F[ABSposChi1])	## removing the diverging element from real part
	Chin_F[ABSposChi2] = 1.0j*sp.imag(Chin_F[ABSposChi2])
	Chia_F[ABSposChi1] = 1.0j*sp.imag(Chia_F[ABSposChi1])
	Chia_F[ABSposChi2] = 1.0j*sp.imag(Chia_F[ABSposChi2])
	return [Chin_F,Chia_F,ABSposChi1,ABSposChi2]


def SelfEnergy(GFn_F,GFa_F,ChiGamma_F,Zn_F):
	""" calculating the dynamical self-energy from Schwinger-Dyson equation,
	    uses fft from scipy.fftpack """
	N = (len(Zn_F)-1)/2
	dz = sp.around(Zn_F[1]-Zn_F[0],8)
	## zero-padding the arrays to double the size
	FD_F = 1.0*sp.concatenate([sp.zeros(3*N+3),sp.ones(N+1)])
	ImGFn_F = sp.concatenate([sp.imag(GFn_F[(len(Zn_F)-1)/2:])\
	,sp.zeros(2*N+3),sp.imag(GFn_F[:(len(Zn_F)-1)/2])])
	ImGFa_F = sp.concatenate([sp.imag(GFa_F[(len(Zn_F)-1)/2:])\
	,sp.zeros(2*N+3),sp.imag(GFa_F[:(len(Zn_F)-1)/2])])
	ImCG_F = sp.concatenate([sp.imag(ChiGamma_F[(len(Zn_F)-1)/2:])\
	,sp.zeros(2*N+3),sp.imag(ChiGamma_F[:(len(Zn_F)-1)/2])])
	## perform convolution/cross-correlation via FFT 
	ftImSEn1_F = sp.conj(fft(FD_F*ImCG_F))*fft(ImGFn_F)*dz
	ftImSEn2_F = -fft(FD_F*ImGFn_F)*sp.conj(fft(ImCG_F))*dz
	ImSEn_F = sp.real(ifft(ftImSEn1_F+ftImSEn2_F))/sp.pi
	ImSEn_F = sp.concatenate([ImSEn_F[3*N+4:],ImSEn_F[:N+1]])
	ftImSEa1_F = sp.conj(fft(FD_F*ImCG_F))*fft(ImGFa_F)*dz
	ftImSEa2_F = -fft(FD_F*ImGFa_F)*sp.conj(fft(ImCG_F))*dz
	ImSEa_F = sp.real(ifft(ftImSEa1_F+ftImSEa2_F))/sp.pi
	ImSEa_F = sp.concatenate([ImSEa_F[3*N+4:],ImSEa_F[:N+1]])
	## find real part from imaginary using KK relations
	Sigman_F = KramersKronigFFT(ImSEn_F) + 1.0j*ImSEn_F
	Sigmaa_F = KramersKronigFFT(ImSEa_F) + 1.0j*ImSEa_F
	return [Sigman_F,Sigmaa_F]


def KramersKronigFFT(ImX_F):
	"""	Hilbert transform used to calculate real part of a function from its imaginary part
		uses piecewise cubic interpolated integral kernel of the Hilbert transform
		assumes that Im X (\infty)=0 
		use only if len(ImX_F)=2**m-1 as it uses fft from scipy.fftpack  """
	N = len(ImX_F)
	A = sp.arange(3,N+1,dtype='float64')	## be careful with the data type
	X1 = 4.0*sp.log(1.5)
	X2 = 10.0*sp.log(4.0/3.0)-6.0*sp.log(1.5)
	## filling the kernel
	Kernel_F = sp.zeros(N-2,dtype='float64')
	Kernel_F = (1-A**2)*((A-2)*sp.arctanh(1.0/(1-2*A))+(A+2)*sp.arctanh(1.0/(1+2*A)))\
	+((A**3-6*A**2+11*A-6)*sp.arctanh(1.0/(3-2*A))+(A+3)*(A**2+3*A+2)*sp.arctanh(1.0/(2*A+3)))/3.0
	Kernel_F = sp.concatenate([-sp.flipud(Kernel_F),sp.array([-X2,-X1,0.0,X1,X2]),Kernel_F])/sp.pi
	## zero-padding the functions for fft
	ImXExt_F = sp.concatenate([ImX_F[(N-1)/2:],sp.zeros(N+2),ImX_F[:(N-1)/2]])
	KernelExt_F = sp.concatenate([Kernel_F[N:],sp.zeros(1),Kernel_F[:N]])
	## performing the fft
	ftReXExt_F = -fft(ImXExt_F)*fft(KernelExt_F)
	ReXExt_F = sp.real(ifft(ftReXExt_F))
	ReX_F = sp.concatenate([ReXExt_F[(3*N+3)/2+1:],ReXExt_F[:(N-1)/2+1]])
	return ReX_F


def KramersKronigFFT_ABS(ImX_F,En_F,ABSpos_F,Res_F):
	"""  Hilbert transform used to calculate real part of a function from its imaginary part
		uses piecewise cubic interpolated integral kernel of the Hilbert transform
		assumes that Im X (\infty)=0 and there are delta-functions at ABSpos_F
		the 1/x polar contributions are added analytically
		use only if len(ImX_F)=2**N-1 as it uses fft from scipy.fftpack  """
	X_F = sp.copy(ImX_F)
	N = len(X_F)
	## removing poles from the function
	if len(ABSpos_F) > 0: 
		for i in range(len(ABSpos_F)):
			X_F[int(ABSpos_F[i])] = 0.0
	A = sp.arange(3,N+1,dtype='float64')	## be careful with the data type
	X1 = 4.0*sp.log(1.5)
	X2 = 10.0*sp.log(4.0/3.0)-6.0*sp.log(1.5)
	## filling the kernel	
	Kernel_F = sp.zeros(N-2,dtype='float64')
	Kernel_F = (1-A**2)*((A-2)*sp.arctanh(1.0/(1-2*A))+(A+2)*sp.arctanh(1.0/(1+2*A)))\
	+((A**3-6*A**2+11*A-6)*sp.arctanh(1.0/(3-2*A))+(A+3)*(A**2+3*A+2)*sp.arctanh(1.0/(2*A+3)))/3.0
	Kernel_F = sp.concatenate([-sp.flipud(Kernel_F),sp.array([-X2,-X1,0.0,X1,X2]),Kernel_F])/sp.pi
	## zero-padding the functions for fft
	ImXExt_F = sp.concatenate([X_F[(N-1)/2:],sp.zeros(N+2),X_F[:(N-1)/2]])
	KernelExt_F = sp.concatenate([Kernel_F[N:],sp.zeros(1),Kernel_F[:N]])
	## performing the fft
	ftReXExt_F = -fft(ImXExt_F)*fft(KernelExt_F)
	ReXExt_F = sp.real(ifft(ftReXExt_F))
	ReX_F = sp.concatenate([ReXExt_F[(3*N+3)/2+1:],ReXExt_F[:(N-1)/2+1]])
	## adding back the 1/x part from the ABS
	if len(ABSpos_F) > 0: 
		if len(ABSpos_F) != 2:
			print('# KramersKronigFFT_ABS: more than two ABS states, exit.')
			exit(-1)
		for i in range(len(ABSpos_F)):
			if ABSpos_F[i] != -1: 
				OneOverX = 1.0/(En_F-En_F[int(ABSpos_F[i])]+1e-10)
				ReX_F = ReX_F+Res_F[i]*OneOverX
	return ReX_F


#####################################################################
# calculation of the interacting Green's function ###################

def GreensFunction(params_F,n,mu,SEn_F,SEa_F,SEnStar_F,SEaStar_F,En_F,A):
	"""	general form the Nambu Green's function and its determinant
	    ABS in gap are not included, must be calculated separately
	    A='band' calculates determinant for band
	    A='gap' calculates determinant for gap """
	[U,Delta,GammaR,GammaL,GammaN,Phi,eps,h] = params_F
	hfe = eps+U*(n-0.5)
	if A == 'band':
		S_F = SFunctionBand(GammaR,GammaL,Delta,En_F)
		D_F = DeltaFunctionBand(GammaR,GammaL,Delta,Phi,En_F)
	else:
		S_F = SFunctionGap(GammaR,GammaL,Delta,En_F)
		D_F = DeltaFunctionGap(GammaR,GammaL,Delta,Phi,En_F)
	Detn_F = (En_F*(1.0+S_F)-hfe-SEn_F)*(En_F*(1.0+S_F)+hfe-SEnStar_F)
	Deta_F = (D_F-U*mu-SEa_F)*(D_F-U*mu-SEaStar_F)
	Det_F = Detn_F-Deta_F
	GFn_F =  (En_F*(1.0+S_F)+hfe-SEnStar_F)/Det_F
	GFa_F = -(D_F-U*mu-SEa_F)/Det_F
	return [GFn_F,GFa_F,Det_F]


def FindABS(Det_F,En_F,Delta):
	"""	determines the energies of ABS as zeroes of GF determinant """
	from subprocess import call
	dE = sp.around(En_F[1]-En_F[0],8)
	dE_dec = int(-sp.log10(dE))
	[EdgePos1,EdgePos2] = FindEdges(En_F,Delta)
	DetG = InterpolatedUnivariateSpline(En_F[EdgePos1+1:EdgePos2],sp.real(Det_F[:]))
	RootsG_F = DetG.roots()
	if len(RootsG_F) == 0:	
		## ABS states too close to gap edges
		## this also happens when using brentq to calculate densities: 
		## it starts from wrong first guess (lower end of bracketing interval)
		ABS_F = sp.array([-Delta+2.0*dE,Delta-2.0*dE])
		ABSpos_F = sp.array([EdgePos1+1,EdgePos2-1])
		Diff_F = sp.array([DetG.derivatives(ABS_F[0])[1],DetG.derivatives(ABS_F[1])[1]])
	elif len(RootsG_F) == 1: 
		## ABS too close to each other?
		print("# Warning: FindABS: only one ABS found: "+str(RootsG_F[0])+", using mirroring to get the other.")
		ABS_F = sp.zeros(2)
		ABSpos_F = sp.zeros(2)
		Diff_F = sp.zeros(2)
		ABS_F = [-sp.fabs(RootsG_F[0]),sp.fabs(RootsG_F[0])]
		for i in range(2):
			ABSpos_F[i] = sp.nonzero(En_F == sp.around(ABS_F[i],dE_dec))[0][0]
			Diff_F[i] = DetG.derivatives(ABS_F[i])[1]
	elif len(RootsG_F) == 2:	
		## two ABS states, ideal case
		ABS_F = sp.zeros(2)
		ABSpos_F = sp.zeros(2)
		Diff_F = sp.zeros(2)
		for i in range(2):
			ABS_F[i] = RootsG_F[i]
			ABSpos_F[i] = sp.nonzero(En_F == sp.around(RootsG_F[i],dE_dec))[0][0]
			Diff_F[i] = DetG.derivatives(RootsG_F[i])[1]
	elif len(RootsG_F) == 6:	
		## four ABS states, other two zeros are in fact poles
		ABS_F = sp.zeros(4)
		ABSpos_F = sp.zeros(4)
		Diff_F = sp.zeros(4)
		D = {0:1, 1:2, 2:3, 3:4}	## maps ABS to zeros of determinant, excludes poles
		for i in range(4):
			ABS_F[i] = RootsG_F[D[i]]
			ABSpos_F[i] = sp.nonzero(En_F == sp.around(RootsG_F[D[i]],dE_dec))[0][0]
			Diff_F[i] = DetG.derivatives(RootsG_F[D[i]])[1]
		## taking two inner poles only
		ABS_F = sp.array([ABS_F[1],ABS_F[2]])
		Diff_F = sp.array([Diff_F[1],Diff_F[2]])
		ABSpos_F = sp.array([int(ABSpos_F[1]),int(ABSpos_F[2])])
	else:
		## this happens when poles are present in DetG and they are fitted in wrong way 
		## taking two inner poles as ABS
		NABS = len(RootsG_F)
		print("# Warning: FindABS: determinant of GF has "+str(NABS)+" poles.")
		ABS_F = sp.zeros(2)
		ABSpos_F = sp.zeros(2)
		Diff_F = sp.zeros(2)
		ABS_F = [-sp.fabs(RootsG_F[NABS/2]),sp.fabs(RootsG_F[NABS/2])]
		for i in range(2):
			ABSpos_F[i] = sp.nonzero(En_F == sp.around(ABS_F[i],dE_dec))[0][0]
			Diff_F[i] = DetG.derivatives(ABS_F[i])[1]
	if len(ABS_F) != 2: 
		## just in case we miss the fact we have more poles than two
		print('# Warning: FindABS: More than two poles???')
	return [ABS_F,Diff_F,ABSpos_F]


def FillGreensFunction(params_F,n,mu,SEn_F,SEa_F,En_F):
	""" calculating the interacting Green's function using the Dyson equation
	    weights of ABS are calculated numerically from determinant, 
	    real parts are recalculated using KK relations """
	[U,Delta,GammaR,GammaL,GammaN,Phi,eps,h] = params_F
	hfe = eps+U*(n-0.5)
	dE = sp.around(En_F[1]-En_F[0],8)
	dE_dec = int(-sp.log10(dE))
	## find special points
	[EdgePos1,EdgePos2] = FindEdges(En_F,Delta)
	SEnStar_F = -sp.flipud(sp.conj(SEn_F))	# hole self-energies
	SEaStar_F =  sp.flipud(sp.conj(SEa_F))
	## calculate GF separately in band and in gap regions
	[GFn1_F,GFa1_F,Det1_F] = GreensFunction(params_F,n,mu,SEn_F[:EdgePos1],SEa_F[:EdgePos1],\
		SEnStar_F[:EdgePos1],SEaStar_F[:EdgePos1],En_F[:EdgePos1],'band')
	[GFn2_F,GFa2_F,Det2_F] = GreensFunction(params_F,n,mu,SEn_F[EdgePos1+1:EdgePos2],SEa_F[EdgePos1+1:EdgePos2],\
		SEnStar_F[EdgePos1+1:EdgePos2],SEaStar_F[EdgePos1+1:EdgePos2],En_F[EdgePos1+1:EdgePos2],'gap')
	[GFn3_F,GFa3_F,Det3_F] = GreensFunction(params_F,n,mu,SEn_F[EdgePos2+1:],SEa_F[EdgePos2+1:],\
		SEnStar_F[EdgePos2+1:],SEaStar_F[EdgePos2+1:],En_F[EdgePos2+1:],'band')
	Det_F = sp.concatenate([Det1_F,sp.zeros(1),Det2_F,sp.zeros(1),Det3_F])
	SigmanStar = UnivariateSpline(En_F[EdgePos1+1:EdgePos2],sp.real(SEnStar_F[EdgePos1+1:EdgePos2]))
	Sigmaa     = UnivariateSpline(En_F[EdgePos1+1:EdgePos2],sp.real(SEa_F[EdgePos1+1:EdgePos2]))
	GFn_F = sp.concatenate([GFn1_F,sp.zeros(1),GFn2_F,sp.zeros(1),GFn3_F])
	GFa_F = sp.concatenate([GFa1_F,sp.zeros(1),GFa2_F,sp.zeros(1),GFa3_F])
	[ABS_F,Diff_F,ABSpos_F] = FindABS(Det2_F,En_F,Delta)
	## the code is being prepared for more than two ABS, but not finished, so we stop here
	if len(ABS_F) !=2: 
		print('# FillGreensFunction: more than two ABS states, exit.')
		exit(-1)
	## calculate the residues
	NABS = len(ABS_F)	# number of ABS
	Sf_F = sp.zeros(NABS,dtype=complex)
	Df_F = sp.zeros(NABS,dtype=complex)
	ResGn_F = sp.zeros(NABS)
	ResGa_F = sp.zeros(NABS)
	for i in range(NABS):
		Sf_F[i] = SFunctionGap(GammaR,GammaL,Delta,ABS_F[i])
		Df_F[i] = DeltaFunctionGap(GammaR,GammaL,Delta,Phi,ABS_F[i])
		ResGn_F[i] =  sp.real((ABS_F[i]*(1.0+Sf_F[i])+hfe-SigmanStar(ABS_F[i]))/Diff_F[i])
		ResGa_F[i] = -sp.real((Df_F[i]-U*mu-Sigmaa(ABS_F[i]))/Diff_F[i])
		GFn_F[int(ABSpos_F[i])] = -1.0j*ResGn_F[i]*sp.pi/dE
		GFa_F[int(ABSpos_F[i])] = -1.0j*ResGa_F[i]*sp.pi/dE
	Res_F = sp.concatenate([ResGn_F,ResGa_F])
	## find real part from imaginary using KK relations
	GFn_F = KramersKronigFFT_ABS(sp.imag(GFn_F),En_F,ABSpos_F,ResGn_F)+1.0j*sp.imag(GFn_F)
	GFa_F = KramersKronigFFT_ABS(sp.imag(GFa_F),En_F,ABSpos_F,ResGa_F)+1.0j*sp.imag(GFa_F)
	## add residues to ABS frequencies
	for i in range(NABS):
		GFn_F[int(ABSpos_F[i])] = -1.0j*ResGn_F[i]*sp.pi/dE
		GFa_F[int(ABSpos_F[i])] = -1.0j*ResGa_F[i]*sp.pi/dE
	return [GFn_F,GFa_F,Det_F,ABS_F,ABSpos_F,Res_F]


def IntDOS(GFn_F,En_F):
	""" the integral over the normal DOS, should be 1.0	"""
	TailL =  sp.imag(GFn_F)[0]*En_F[0]/sp.pi	# left tail
	TailR = -sp.imag(GFn_F)[-1]*En_F[-1]/sp.pi	# right tail
	return -trapz(sp.imag(GFn_F),En_F)/sp.pi + TailL + TailR


def MSumsInt(params_F,n,mu,SEn_F,SEa_F,Zn_F):
	"""	calculating Matsubara sums used in calculating n and mu from interacting GF
	    returning three sums M, n = M[1]/(1-U*M[0]), mu = -M[2]/(1-U*M[0])
	    this approach is numerically more precise than integrating the GF """
	[U,Delta,GammaR,GammaL,GammaN,Phi,eps,h] = params_F
	ed = eps-U/2.0
	dz = sp.around(Zn_F[1] - Zn_F[0],8)
	[EdgePos1,EdgePos2] = FindEdges(Zn_F,Delta)
	SEnStar_F = -sp.flipud(sp.conj(SEn_F))	# hole self-energies
	SEaStar_F =  sp.flipud(sp.conj(SEa_F))
	[GFn1_F,GFa1_F,Det1_F] = GreensFunction(params_F,n,mu,SEn_F[:EdgePos1],SEa_F[:EdgePos1],\
		SEnStar_F[:EdgePos1],SEaStar_F[:EdgePos1],Zn_F[:EdgePos1],'band')
	[GFn2_F,GFa2_F,Det2_F] = GreensFunction(params_F,n,mu,SEn_F[EdgePos1+1:EdgePos2],SEa_F[EdgePos1+1:EdgePos2],\
		SEnStar_F[EdgePos1+1:EdgePos2],SEaStar_F[EdgePos1+1:EdgePos2],Zn_F[EdgePos1+1:EdgePos2],'gap')
	[GFn3_F,GFa3_F,Det3_F] = GreensFunction(params_F,n,mu,SEn_F[EdgePos2+1:],SEa_F[EdgePos2+1:],\
		SEnStar_F[EdgePos2+1:],SEaStar_F[EdgePos2+1:],Zn_F[EdgePos2+1:],'band')
	[ABS_F,Diff_F,ABSpos_F] = FindABS(Det2_F,Zn_F,Delta)
	Band_F = Zn_F[:EdgePos1]
	if len(Band_F) > 0: ## Delta smaller than energy minimum
		S_F = SFunctionBand(GammaR,GammaL,Delta,Band_F)
		D_F = DeltaFunctionBand(GammaR,GammaL,Delta,Phi,Band_F)
		Nom2_F = Band_F*(1.0+S_F)+ed-SEnStar_F[:EdgePos1]
		Nom3_F = D_F-SEa_F[:EdgePos1]
		Int1_F = sp.imag(1.0/Det1_F)
		Int2_F = sp.imag(Nom2_F/Det1_F)
		Int3_F = sp.imag(Nom3_F/Det1_F)
		Tail1 = -Int1_F[0]*Zn_F[0]/2.0		## behaves as 1/x^3
		Tail2 = -Int2_F[0]*Zn_F[0]			## behaves as 1/x^2
		Tail3 = -Int3_F[0]*Zn_F[0]/2.0		## behaves as 1/x^3
		Head1 = 0.5*dz*Int1_F[-1]
		Head2 = 0.5*dz*Int2_F[-1]
		Head3 = 0.5*dz*Int3_F[-1]
	else: ## large Delta
		Int1_F = Int2_F = Int3_F = sp.zeros(1)
		Tail1 = Tail2 = Tail3 = 0.0
		Head1 = Head2 = Head3 = 0.0
	if len(ABS_F) == 2:
		Swzero = SFunctionGap(GammaR,GammaL,Delta,ABS_F[0])
		Dwzero = DeltaFunctionGap(GammaR,GammaL,Delta,Phi,ABS_F[0])
		Res2 = ABS_F[0]*(1.0+Swzero)+ed-SEnStar_F[int(ABSpos_F[0])]
		Res3 = Dwzero-SEa_F[int(ABSpos_F[0])]
		MSum1R = -(trapz(Int1_F,Band_F)+Head1+Tail1)/sp.pi+1.0/Diff_F[0]
		MSum2R = -(trapz(Int2_F,Band_F)+Head2+Tail2)/sp.pi+Res2/Diff_F[0]
		MSum3R = -(trapz(Int3_F,Band_F)+Head3+Tail3)/sp.pi+Res3/Diff_F[0]
	else:
		print('# MSumsInt: more than two ABS states, exit')
		exit(-1)
	'''
	elif len(ABS_F) == 4:  ## not yet implemented
		print '# Warning: Four ABS states, taking the two inner only!'
		Swzero1 = SFunctionGap(GammaR,GammaL,Delta,ABS_F[0])
		Swzero2 = SFunctionGap(GammaR,GammaL,Delta,ABS_F[1])
		Dwzero1 = DeltaFunctionGap(GammaR,GammaL,Delta,Phi,ABS_F[0])
		Dwzero2 = DeltaFunctionGap(GammaR,GammaL,Delta,Phi,ABS_F[1])
		Res21 = ABS_F[0]*(1.0+Swzero1)+ed-SEnStar_F[int(ABSpos_F[0])]
		Res22 = ABS_F[1]*(1.0+Swzero2)+ed-SEnStar_F[int(ABSpos_F[1])]
		Res31 = Dwzero1-SEa_F[(ABSpos_F[0])]
		Res32 = Dwzero2-SEa_F[(ABSpos_F[1])]
		MSum1R = -(trapz(Int1_F,Band_F)+Head1+Tail1)/sp.pi+1.0/Diff_F[0]+1.0/Diff_F[1]
		MSum2R = -(trapz(Int2_F,Band_F)+Head2+Tail2)/sp.pi+Res21/Diff_F[0]+Res22/Diff_F[1]
		MSum3R = -(trapz(Int3_F,Band_F)+Head3+Tail3)/sp.pi+Res31/Diff_F[0]+Res32/Diff_F[1]
	'''
	return sp.real_if_close([MSum1R,MSum2R,MSum3R])


def ElectronDensity(params_F,n,mu,SEn_F,SEa_F,En_F):
	"""	calculating n from Matsbara sums MSumsInt """
	U = params_F[0]
	MSums_F = MSumsInt(params_F,n,mu,SEn_F,SEa_F,En_F)
	n = sp.real_if_close(MSums_F[1]/(1.0 - U*MSums_F[0]))
	if sp.imag(n) > 0.0: print('# Warning: non-zero imag. part of n: {0: .8f}'.format(float(sp.imag(n))))
	return sp.real(n)


def CooperPairDensity(params_F,n,mu,SEn_F,SEa_F,En_F):
	""" calculating mu from Matsbara sums MSumsInt """
	U = params_F[0]
	MSums_F = MSumsInt(params_F,n,mu,SEn_F,SEa_F,En_F)
	mu =  sp.real_if_close(-MSums_F[2]/(1.0 - U*MSums_F[0]))
	if sp.imag(mu) > 0.0: print('# Warning: non-zero imag. part of mu: {0: .8f}'.format(float(sp.imag(mu))))
	return sp.real(mu)


def JosephsonCurrent(params_F,GFa_F,En_F,ResGa,wzero):
	""" calculates the Josephson current, separated into the band part and the gap (ABS) part """
	[U,Delta,GammaR,GammaL,GammaN,Phi,eps,h] = params_F
	[EdgePos1,EdgePos2] = FindEdges(En_F,Delta)
	PreFac   = Delta*(GammaL+GammaR)*sp.sin(Phi/2.0)
	Int_F    = sp.real(GFa_F[:EdgePos1])/(sp.sqrt(En_F[:EdgePos1]**2-Delta**2))
	if Delta < sp.fabs(En_F[0]):	## finite gap
		BandPart = PreFac*simps(Int_F,En_F[:EdgePos1])/sp.pi
	else: BandPart = 0.0
	GapPart  = PreFac*ResGa/sp.sqrt(Delta**2-wzero**2)
	return sp.array([BandPart,GapPart])


def SmoothWings(En_F,X_F,Emax):
	""" smooths a function on given energy interval, not used anymore """
	de = sp.around(En_F[1]-En_F[0],8)
	de_dec = int(-sp.log10(de))
	PosX1 = sp.nonzero(En_F == sp.around(-Emax,de_dec))[0][0]
	PosX2 = sp.nonzero(En_F == sp.around(Emax,de_dec))[0][0]
	print('# SmoothWings: start...')
	Rex1 = UnivariateSpline(En_F[:PosX1],sp.real(X_F[:PosX1]),s=1e-3)
	Rex2 = UnivariateSpline(En_F[PosX2:],sp.real(X_F[PosX2:]),s=1e-3)
	print('# SmoothWings: end.')
	X_F = sp.concatenate([Rex1(En_F[:PosX1]),sp.real(X_F[PosX1:PosX2]),Rex2(En_F[PosX2:])])+1.0j*sp.imag(X_F)
	return X_F


def FitTail(En_F,X_F,Emin,Emax,parity):
	""" fits tails of real part of a function on -Emax:-Emin and Emin:Emax
    	uses first two orders of either odd or even symmetry """
	from scipy.optimize import curve_fit
	de = sp.around(En_F[1]-En_F[0],8)
	de_dec = int(-sp.log10(de))
	Pos1 = sp.nonzero(En_F == sp.around(-Emax,de_dec))[0][0]
	Pos2 = sp.nonzero(En_F == sp.around(-Emin,de_dec))[0][0]
	Pos3 = sp.nonzero(En_F == sp.around( Emin,de_dec))[0][0]
	Pos4 = sp.nonzero(En_F == sp.around( Emax,de_dec))[0][0]
	if parity == 'odd':
		f = lambda x,a,b: a/x + b/x**3
	elif parity == 'even':
		f = lambda x,a,b: a/x**2 + b/x**4
	try:
		params1_F = curve_fit(f,En_F[Pos1:Pos2],sp.real(X_F[Pos1:Pos2]))[0]
		params2_F = curve_fit(f,En_F[Pos3:Pos4],sp.real(X_F[Pos3:Pos4]))[0]
		X_F = sp.concatenate([f(En_F[:Pos1],params1_F[0],params1_F[1]),\
		sp.real(X_F[Pos1:Pos3]),f(En_F[Pos3:],params2_F[0],params2_F[1])])+1.0j*sp.imag(X_F)
	except RuntimeError:
		print('# FitTail: Warning: tails cannot be fitted.')
	return X_F

## squadlib2.py end ##

