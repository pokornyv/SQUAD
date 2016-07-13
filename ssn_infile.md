Description of the *ssn.in* file
================================

The file *ssn.in* defines various parameters used by the
*ssn_second.py* code that calculates properties of a system with two 
sc and one normal electrodes.  
  
##Description:##  

### [params] section: ###

- M : the energy axis contains 2^M+1 points. Default: 20  
- dE : discretization of the energy axis. Default: 1e-4  
- rootf : method used in calcualtion of the densities n=<d+d> and mu=<d+d+>. Options are *brentq* and *fixed_point*.  
- ConvN : convergence criterium used in self-consitent calculation of n and mu in 2nd-order PT solution. Default: 1e-4  
- ConvX : convergence criterium used in calculation n and mu using method defined in *rootf*. Default: 1e-5  
- ConvHF : convergence criterium used in calculation of n and mu in the Hartree-Fock solution          :  1e-6  
- MuMin : minimum of the separated interval where we search for mu. Used if *rootf=brentq*. Default: -2.0  
- MuMax : maximum of the separated interval where we search for mu. Used if *rootf=brentq*. Default:  2.0  
- HF_max_iter : maximum number of iterations in the Hartree-Fock solver. Default: 10000
- offset_x : offset of the energies used to avoid poles in functions (e.g. gap edge). Default: 1e-12

### [IO] section: ###
  
- WriteIO :  0/1 switch whether we want to write calculation details to standard output (e.g. screen)
- WriteFile_HFGF   :  0/1 switch whether the file with HF Green function should be written  
- WriteFile_Bub    :  0/1 switch whether the file with HF bubble should be written  
- WriteFile_2ndSE  :  0/1 switch whether the file with 2nd-order self-energy should be written  
- WriteFile_2ndGF  :  0/1 switch whether the file with 2nd-order Green function should be written  
  
- EmaxFiles : maximum of the energy window for file output. Default: 10  
- EstepFiles : step in energies for file output. The values will be written with EstepFiles x dE step. Default: 10  