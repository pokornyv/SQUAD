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
- Modify parameters in *squad.in*. The default values are OK for most cases, change only if needed.  
- Run `python secondPT.py <float U> <float Δ> <float ΓR> <float ΓLR> <float ε> <float P>`  
  where *python* is the alias for the python 2 interpreter, ΓLR is GammaL/GammaR and P=Φ/π.  
- Check standard output for solution. Parameters like ABS energy, densities, ABS residues and the supercurrent are printed.  
- The Green function and the self-energy is printed to files (check *\[IO\]* section of *squad.in* for details).  


#### Basic features:
- based on the diagrammatic perturbation theory in Coulomb interaction U in Nambu formalism  
- calculates properties within the second order perturbation theory (2ndPT) at zero temperature  
- weak-coupling expansion, works well for small U/Γ  
- real-frequency representation (exact analytic continuation from Matsubara formalism)  
- direct access to spectral functions and Andreev bound state (ABS) energies  
- reliable description of the 0-phase and the 0-π phase boundary  
- π-phase is not accesible due to the double-degenerate ground state  

#### List of files:
- *squadlib1.py* - library of general functions and the Hartree-Fock solver  
- *squadlib2.py* - library of functions for calculating 2ndPT  
- *secondPT.py* - main code to calculate 2nd order PT results for a system with two sc electrodes  
- *params.py* - script to read ss.in file  
- *squad.in* - parameter file for *second_sc.py*, read by *params.py*, described in *infile.md*  

#### Known issues:
- Large number of points on the energy axis can lead to `RuntimeWarning: invalid value encountered in power`
warning. We need large arrays for the Hilbert transform, where it can 
hit the limit of 2**63 of signed int. Reduce the number of points (*M* in *squad.in*).  
- If the energy of the Andreev bound state is close to the gap edge Δ, 'AndreevEnergy()' solver
can fail. Please change the *ABSinit_val* parameter in *squad.in*.  
- Close to the 0-π transition, the solver can become unstable as it cannot describe the π-phase.
Check if the densities n and mu are real numbers and pay attention to warnings. Also, the residues
may behave erratically.  
- Beyond the 0-π transition, the solver gives zero ABS energies or 
it can end up in an unphysical "false-π" phase. It is marked
by the change of the sign of the residues in anomalous Green function and (therefore) negative mu. 
Pay attention to warnings.  

#### TODO list:
- [ ] include infinite-Delta solution  
- [ ] include magnetic field to study the splitting of the Andreev states  

#### few solutions as a benchmark:
```
  U     GammaR  GammaL  eps     Phi/pi  wABS        n           mu          ResGn1      ResGn2      ResA1

- eps - dependence:
 1.000	 0.500	 0.500	-0.500	 0.500	 0.30510	 0.66949	 0.17922	 0.38401	 0.12737	 0.22118
 1.000	 0.500	 0.500	 0.000	 0.500	 0.23882	 0.50000	 0.21346	 0.25985	 0.25985	 0.25991
 1.000	 0.500	 0.500	 0.500	 0.500	 0.30511	 0.33054	 0.17922	 0.12737	 0.38402	 0.22118

- phi-dependence:
 2.000	 0.500	 0.500	 0.000	 0.000	 0.26748	 0.50000	 0.18751	 0.26556	 0.26556	 0.26568
 2.000	 0.500	 0.500	 0.000	 0.250	 0.22462	 0.50000	 0.18726	 0.26011	 0.26011	 0.26034
 2.000	 0.500	 0.500	 0.000	 0.500	 0.10833	 0.50000	 0.18640	 0.24832	 0.24832	 0.24966
```
Results may vary slightly, depending on the parameters in *squad.in*.  

#### References:
1. M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Phys. Rev. B 93, 024523 (2016).  
2. M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Sci. Rep. 5, 8821 (2015).  

