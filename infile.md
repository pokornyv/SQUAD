Description of the *squad.in* file
==================================

File *squad.in* defines parameters used by *secondPT.py* code.  
  
## Description

### [params] section
  
- M : the energy axis contains 2^M+1 points. Default: 21  
- dE : discretization of the energy axis. Default: 1e-4  
- rootf : method for calculation of the densities n=\<d\+d\> and μ=\<d\+d\+\>. Options are *brentq* and *fixed_point*.  
- ConvN : convergence criterion for the self-consitent calculation of n and μ in 2nd-order PT solution. Default: 1e-4  
- ConvX : convergence criterion for the calculation n and μ using method defined in *rootf*. Default: 1e-5  
- ConvHF : convergence criterion for the calculation of n and μ in the Hartree-Fock solution. Default: 1e-6  
- MuMin : minimum of the separated interval where we search for μ. Used if *rootf=brentq*. Default: -2  
- MuMax : maximum of the separated interval where we search for μ. Used if *rootf=brentq*. Default:  2  
- ABSinit_val : initial value to start a fixed-point calculation of the ABS energy is *ABSinit_val x Delta*. Default: 0.99  
- HF_max_iter : maximum number of iterations for the Hartree-Fock solver. Default: 10000  
- offset_x : offset of the energies used to avoid poles in functions (e.g. gap edges). Default: 1e-12  

### [IO] section

0/1 (yes/no) switches:
- WriteIO          :  write calculation details to standard output  
- WriteFile_HFGF   :  write *HF_green.dat* file with HF Green function  
- WriteFile_Bub    :  write *HF_bubbles.dat* file with HF two-particle bubble  
- WriteFile_2ndSE  :  write *2nd_SE.dat* file with 2nd-order self-energy  
- WriteFile_2ndGF  :  write *2nd_green.dat* file with 2nd-order Green function  

Other output parameters:
- EmaxFiles : maximum of the energy window for output. Default: 20.0  
- EstepFiles : energy step for output. Values will be written with (EstepFiles x dE) step. Default: 10  

