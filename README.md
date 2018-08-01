SQUAD
=====
#### Description:
**SQUAD** (**su**perconducting **qua**ntum **d**ot) is a python2 code to calculate spectral 
and transport properties of a single-level quantum dot connected to two BCS superconducting (sc) leads.
  
SQUAD requires SciPy libraries (www.scipy.org) and it was developed and tested using python 2.7.5 and SciPy 0.12 and 0.19.  
  
These codes are subject to changes and require some optimization and cleanup. Please use with caution and at own risk.
Please check Ref. 1 and 2 before start using the code. In a case you find a bug or need help using the code, 
please create an issue on GitHub (preffered) or contact me at pokornyv@fzu.cz. Please always include the 
commit hash prefix (e.g. *7719c9a*) in all communication. 

SQUAD is a free software distributed under the GPL license.

#### Project homepage:
https://github.com/pokornyv/SQUAD

#### Basic features:
- based on the diagrammatic perturbation theory in Coulomb interaction U in Nambu formalism  
- calculates properties within the second order perturbation theory (2ndPT) at zero temperature  
- weak-coupling expansion, works well for small U/Γ  
- real-frequency representation (exact analytic continuation from Matsubara formalism)  
- direct access to spectral functions and Andreev bound state (ABS) energies  
- reliable description of the 0-phase and the 0-π phase boundary  
- π-phase is not accesible due to the double-degenerate ground state  

#### Usage:
- Modify parameters in *squad.in*, if needed. The default values are sufficient for most cases.  
- Run `python secondPT.py <float U> <float Δ> <float ΓR> <float ΓLR> <float ε> <float P>`  
where *python* is the alias for the python 2 interpreter, ΓLR is GammaL/GammaR, 
ε is the local energy (gate voltage) w.r.t. half-filling (so ε=0 represents a half-filled dot) and P=Φ/π.  
- Check standard output for the solution. Parameters like ABS energy, densities n=\<d\+d\> and μ=\<d\+d\+\>, 
ABS residues and the supercurrent are printed.  
- The Green function and the self-energy can be printed to files (check \[IO\] section of *squad.in* for details).  

#### List of files:
- *squadlib1.py* - library of general functions and the Hartree-Fock solver  
- *squadlib2.py* - library of functions for calculating 2ndPT  
- *secondPT.py* - main code to calculate 2nd order PT results for a system with two sc electrodes  
- *params.py* - script to read *squad.in* file  
- *squad.in* - parameter file for *secondPT.py*, read by *params.py*, described in *infile.md*  
- *infile.md* - description of the *squad.in* file  
- *README.md* - this file  

#### Known issues:
- Large number of points on the energy axis can lead to `RuntimeWarning: invalid value encountered in power`
warning. We need large arrays for the Hilbert transform, where it can 
hit the limit of 2**63 of signed int. Reduce the number of points (*M* in *squad.in*).
For most applications, *M=20* with *dE=1e-4* is sufficient to get results correct up to 3 decimal places,
except the vicinity of the phase transition, where ABS energy can fall below dE.  
- If the energy of the Andreev bound state is close to the gap edge Δ, `AndreevEnergy()` solver
can fail. Please change the *ABSinit_val* parameter in *squad.in* to some reasonable value.  
- Close to the 0-π transition, the solver can become unstable as it cannot describe the π-phase.
Check if the densities n and μ are real numbers and pay attention to warnings. Also, the residues
may behave erratically, having negative effect on the Josephson current and the densities n and μ.  
- Beyond the 0-π transition, the solver gives zero ABS energies or it can end up in an unphysical "false-π" phase. It is marked
by the change of the sign of the residues in anomalous Green function and (therefore) negative μ. 
Pay attention to warnings.  
- There is an instability in the calculation for Φ=π (P=1). Please use e.g. P=0.999 until this issue is solved.

#### TODO list:
- [ ] Fix the instability at Φ=π.  
- [ ] Fix the instability in calculating residues close to the quantum phase transition.
- [ ] Include infinite-Delta solution. This requires modifications to the Hartree-Fock solver.  
- [ ] Include magnetic field to study the splitting of the Andreev bound states.  

#### few solutions as a benchmark:
```
  U     Delta   GammaR  GammaL  eps     Phi/pi  wABS        n           mu          ResGn1      ResGn2      ResGa1      JC

- eps-dependence for U = 1:
 1.000	 1.000	 0.500	 0.500	-1.000	 0.500	 0.45552	 0.76912	 0.12708	 0.40992	 0.06238	 0.15991	 0.08007
 1.000	 1.000	 0.500	 0.500	-0.500	 0.500	 0.30510	 0.66949	 0.17922	 0.38401	 0.12737	 0.22118	 0.11980
 1.000	 1.000	 0.500	 0.500	 0.000	 0.500	 0.23882	 0.50000	 0.21346	 0.25985	 0.25985	 0.25991	 0.14665
 1.000	 1.000	 0.500	 0.500	 0.500	 0.500	 0.30510	 0.33053	 0.17922	 0.12737	 0.38402	 0.22118	 0.11980
 1.000	 1.000	 0.500	 0.500	 1.000	 0.500	 0.45553	 0.23091	 0.12708	 0.06238	 0.40992	 0.15990	 0.08007

- phi-dependence for U = 2:
 2.000	 1.000	 0.500	 0.500	 0.000	 0.000	 0.26748	 0.50000	 0.18751	 0.26556	 0.26556	 0.26568	 0.00000
 2.000	 1.000	 0.500	 0.500	 0.000	 0.100	 0.26048	 0.50000	 0.18744	 0.26459	 0.26459	 0.26473	 0.03046
 2.000	 1.000	 0.500	 0.500	 0.000	 0.200	 0.23982	 0.50000	 0.18731	 0.26191	 0.26191	 0.26210	 0.06045
 2.000	 1.000	 0.500	 0.500	 0.000	 0.300	 0.20659	 0.50000	 0.18706	 0.25793	 0.25793	 0.25823	 0.08988
 2.000	 1.000	 0.500	 0.500	 0.000	 0.400	 0.16218	 0.50000	 0.18678	 0.25330	 0.25330	 0.25387	 0.11852
 2.000	 1.000	 0.500	 0.500	 0.000	 0.500	 0.10833	 0.50000	 0.18640	 0.24832	 0.24832	 0.24966	 0.14629
 2.000	 1.000	 0.500	 0.500	 0.000	 0.600	 0.04669	 0.50000	 0.18593	 0.24321	 0.24321	 0.24839	 0.17494
 2.000	 1.000	 0.500	 0.500	 0.000	 0.650	 0.01345	 0.50000	 0.18566	 0.24056	 0.24056	 0.26463	 0.20312
 2.000	 1.000	 0.500	 0.500	 0.000	 0.660	 0.00661	 0.50000	 0.18561	 0.24004	 0.24004	 0.29224	 0.22983
 2.000	 1.000	 0.500	 0.500	 0.000	 0.670	 0.00000	 0.50000	 0.18529	 0.23949	 0.23949	 -------     -------

```
Results may vary slightly, depending on the parameters in *squad.in*.  

#### References:
1. M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Phys. Rev. B 93, 024523 (2016).  
2. M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Sci. Rep. 5, 8821 (2015).  

