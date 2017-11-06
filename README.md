SQUAD
=====
#### Description:
**SQUAD** (**su**perconducting **qua**ntum **d**ot) is a python code to calculate spectral 
and transport properties of a single-level quantum dot connected to two BCS superconducting (sc) leads.
Third, non-sc lead can be also connected (codes for calculations with normal electrode have *ssn* in their names),
but this solution works correctly only for very small interaction strength.
  
SQUAD requires SciPy libraries (www.scipy.org) and it was developed and tested using python 2.7.5 and SciPy 0.12.  
  
These codes are subject to big changes and require some optimization and cleanup. Please use with caution. 
In case of interest please contact me at pokornyv@fzu.cz.

#### Project homepage:
https://github.com/pokornyv/SQUAD

#### Usage:
- `python second_sc.py <float U> <float Delta> <float GammaR> <float > <float GLR> <float eps> <float P>`  
- `python ssn_second.py <float U> <float Delta> <float GammaR> <float > <float GLR> <GammaN> <float eps> <float P>`  

where *python* is the alias for the python 2 interpreter, GLR is GammaL/GammaR and P=Phi/pi.

#### Basic features:
- based on the diagrammatic perturbation theory in Coulomb interaction U in Nambu formalism  
- calculates properties within the second order perturbation theory (PT)  
- weak-coupling expansion, works well for small U/Gamma  
- real-frequency representation (analytic continuation from Matsubara formalism)  
- direct access to spectral functions and Andreev bound states  
- reliable description of the zero-phase and the zero-pi phase boundary  
- pi-phase is not accesible due to the double-degenerate ground state  

#### Solution for infinite Delta:
The code solving the ssn system allows to calculate an approximation of the system with infinitely large gap. Setting Delta large 
(e.g. 1e7 and larger) usually gives a solution numerically correct to 4-5 decimal places. Be sure you set Delta larger
thatn the maximum of the energy axis to get rid of the band contributions.  

#### List of files:
- *squadlib1.py* - library of general functions and the HF solver  
- *squadlib2.py* - library of functions for calculating 2nd order PT  
- *ssnlib.py* - library of functions for calculating ssn results  
- *second_sc.py* - main code to calculate 2nd order PT results for a system with two sc electrodes  
- *ssn_second.py* - main code to calculate 2nd order PT results in a system with one normal and two
sc electrodes  
- *params_ss.py* - script to read ss.in file  
- *params_ssn.py* - script to read ssn.in file  
- *ss.in* - parameter file for *second_sc.py*, read by *params_ss.py*, described in *ssn_infile.md*  
- *ssn.in* - parameter file for *ssn_second.py*, read by *params_ssn.py*, described in *ssn_infile.md*  

#### TODO list:
- [ ] include infinite-Delta solution in the two-terminal (ss) code
- [ ] include magnetic field to study the splitting of the Andreev states
- [ ] calculate reasonable guess for the initial condition for fixed_point in 'AndreevEnergy()'

#### References:
- M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Phys. Rev. B 93, 024523 (2016).  
- M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Sci. Rep. 5, 8821 (2015).  

