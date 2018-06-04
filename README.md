SQUAD
=====
#### Description:
**SQUAD** (**su**perconducting **qua**ntum **d**ot) is a python2 code to calculate spectral 
and transport properties of a single-level quantum dot connected to two BCS superconducting (sc) leads.
  
SQUAD requires SciPy libraries (www.scipy.org) and it was developed and tested using python 2.7.5 and SciPy 0.12.  
  
These codes are subject to changes and require some optimization and cleanup. Please use with caution and at own risk.
Please check Ref. 1 and 2 before start using the code. In a case you find a bug or need help using the code, 
please create an issue on GitHub (preffered) or contact me at pokornyv@fzu.cz. Please always include the 
commit hash prefix (e.g. *8aa17d8*) in all communication. 

SQUAD is free software distributed under the GPL license.

#### Project homepage:
https://github.com/pokornyv/SQUAD

#### Usage:
- `python secondPT.py <float U> <float Δ> <float ΓR> <float ΓLR> <float ε> <float P>`  

where *python* is the alias for the python 2 interpreter, ΓLR is GammaL/GammaR and P=Φ/π.

#### Basic features:
- based on the diagrammatic perturbation theory in Coulomb interaction U in Nambu formalism  
- calculates properties within the second order perturbation theory (2ndPT)  
- weak-coupling expansion, works well for small U/Γ  
- real-frequency representation (exact analytic continuation from Matsubara formalism)  
- direct access to spectral functions and Andreev bound state frequencies  
- reliable description of the 0-phase and the 0-π phase boundary  
- π-phase is not accesible due to the double-degenerate ground state  

#### List of files:
- *squadlib1.py* - library of general functions and the Hartree-Fock solver  
- *squadlib2.py* - library of functions for calculating 2ndPT  
- *secondPT.py* - main code to calculate 2nd order PT results for a system with two sc electrodes  
- *params.py* - script to read ss.in file  
- *squad.in* - parameter file for *second_sc.py*, read by *params.py*, described in *infile.md*  

#### TODO list:
- [ ] include infinite-Delta solution 
- [ ] include magnetic field to study the splitting of the Andreev states
- [ ] calculate reasonable guess for the initial condition for fixed_point in 'AndreevEnergy()'

#### References:
1. M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Phys. Rev. B 93, 024523 (2016).  
2. M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Sci. Rep. 5, 8821 (2015).  

