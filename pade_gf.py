# analytic continuation of Matsubara GF to real axis using Pade approximants
# Vladislav Pokorny; 2015; pokornyv@fzu.cz

import scipy as sp
from pytriqs.gf.local import GfImFreq,GfReFreq,TailGf

filename = 'GF_mats.dat'	# input file with Matsubara GF, columns: n, w_n, Re G, Im G
emax     = 8.0				# real frequency window maximum
de       = 0.01				# step in real frequencies
L        = 200				# number od points from which we do the continuation
eta      = 0.01				# imaginary shift, w = x + i*eta

Data_F = sp.loadtxt(filename)	# load data file to array [NMats:4]

# define and fill the Matsubara GF
NMats  = len(Data_F)	
dwn = sp.fabs(Data_F[1][1]-Data_F[0][1])
beta = 2.0*sp.pi/dwn
print '# beta = {0: .3f}'.format(beta)

GFi = GfImFreq(indices = [0], beta = beta, n_points = NMats)
for i in range(NMats):
	GFi.data[i] = Data_F[i][2]+1.0j*Data_F[i][3]

# fitting the tail of the Matsubara GF, up to 6th order
fixed_tail     = TailGf(1,1,3,-1)
fixed_tail[-1] = sp.zeros([1,1])
fixed_tail[ 0] = sp.zeros([1,1])
fixed_tail[ 1] = sp.eye(1)
GFi.fit_tail(fixed_tail,6,int(0.8*NMats),NMats)

# define and fill the real-time GF
NRealP = int(2.0*emax/de)
GFr    = GfReFreq(indices = [0], window = (-emax,emax), n_points = NRealP)

GFr.set_from_pade(GFi, n_points = L, freq_offset = eta)

# fill the real frequency mesh for output
Freq_F = sp.zeros(NRealP)
k = 0
for w in GFr.mesh:
	Freq_F[k] = sp.real(w)
	k = k+1

# print output
for i in range(len(GFr.mesh)):
	print '{0: .5f}\t{1: .8f}\t{2: .8f}\t'\
	.format(float(Freq_F[i]),float(sp.real(GFr.data[i][0][0])),float(sp.imag(GFr.data[i][0][0])))

