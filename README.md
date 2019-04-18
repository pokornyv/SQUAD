SQUAD
=====
#### Description:
**SQUAD** (**su**perconducting **qua**ntum **d**ot) is a python2 code to calculate spectral 
and transport properties of a single-level quantum dot connected to two BCS superconducting (sc) leads
using second-order perturbation theory (2ndPT) as described in Refs. [1-2].
  
SQUAD requires SciPy [SciPy](https://www.scipy.org) libraries and it was developed and 
tested using python 2.7.5 and SciPy 0.19 and 1.1.
  
Theis codes is subject to constant changes and optimization. Please use at own risk.
Please consult Ref. [1-2] before using the code. In a case you find a bug or need help using the code, 
please create an issue on GitHub (preffered) or contact me at pokornyv@fzu.cz. Please always include the 
commit hash prefix (e.g. *7719c9a*) in all communication.

Please cite Ref. [1] and/or [2] when publishing data calculated using SQUAD.

SQUAD is a free software distributed under the GPL license.

#### Project homepage:
[github.com/pokornyv/SQUAD](https://github.com/pokornyv/SQUAD)

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
- Run `python secondPT.py <U> <Δ> <ΓR> <ΓLR> <ε> <P>`  
where *python* is the alias for the python 2 interpreter, ΓLR is ΓL/ΓR, 
ε is the local energy (gate voltage) w.r.t. half-filling (so ε=0 represents a half-filled dot) and P=Φ/π.  
- Check standard output for the solution. Parameters like ABS energy, densities n=\<d\+d\> and μ=\<d\+d\+\>, 
ABS residues and the supercurrent are printed.  
- The Green function and the self-energy can be printed to files (check \[IO\] section of *squad.in* for details).  

#### List of files:
- *secondPT.py* - main code to calculate 2nd order PT results for a system with two sc electrodes  
- *config_squad.py* - script to read *squad.in* control file and set up global variables  
- *squadlib1.py* - library of general functions and the Hartree-Fock solver  
- *squadlib2.py* - library of functions for calculating 2ndPT  
- *squad.in* - parameter file for *secondPT.py*, described in *infile.md*  
- *infile.md* - description of the *squad.in* file  
- *LICENSE* - a copy of the GNU General Public License  
- *README.md* - this document  

#### Output files:
The output files are written if the corresponding flag is set in *squad.in*. The columns of
the file are always E, Re Xn, Im Xn, Re Xa, Im Xa. The energy interval and density are controlled 
by *EmaxFiles* and *EstepFiles* flags in *squad.in*.

- *HF_green.dat* - Hartree-Fock Green function, set flag *Write_HFGF : 1*  
- *HF_bubbles.dat* - Hartree-Fock bubbles, set flag *Write_Bubble : 1*  
- *2nd_SE.dat* - 2ndPT dynamic self-energy, set flag *Write_2ndSE : 1* Static (HF) parts must be added by hand to real parts.  
- *2nd_green.dat* - 2ndPT Green function, set flag *Write_2ndGF : 1*  

#### Known issues:
- Real parts of bubbles and the 2ndPT Green function are noisy at high energies. This is a result
of the Hilbert transform in Kramers-Kronig relations. This has no effect on the spectral functions. 
Fix is on the way.
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
- [x] Fix the instability in calculating residues close to the quantum phase transition.
- [ ] Fix the instability at Φ=π.  
- [ ] Include magnetic field to study the splitting of the Andreev bound states.  

#### few solutions as a benchmark:
```
  U     Delta   GammaR  GammaL  eps     Phi/pi  wABS        n           mu          ResGn1      ResGn2      ResGa1      JC

- eps-dependence for U = 1:
 1.000	 1.000	 0.500	 0.500	-2.000	 0.5000	 0.75747	 0.85867	 0.06675	 0.25823	 0.02105	 0.07373	 0.03713
 1.000	 1.000	 0.500	 0.500	-1.500	 0.5000	 0.62032	 0.82438	 0.09032	 0.35582	 0.03485	 0.11136	 0.05336
 1.000	 1.000	 0.500	 0.500	-1.000	 0.5000	 0.45551	 0.76912	 0.12708	 0.40992	 0.06238	 0.15991	 0.08006
 1.000	 1.000	 0.500	 0.500	-0.500	 0.5000	 0.30509	 0.66948	 0.17922	 0.38402	 0.12739	 0.22118	 0.11978
 1.000	 1.000	 0.500	 0.500	 0.000	 0.5000	 0.23881	 0.50000	 0.21345	 0.25988	 0.25988	 0.25988	 0.14661
 1.000	 1.000	 0.500	 0.500	 0.500	 0.5000	 0.30509	 0.33054	 0.17922	 0.12738	 0.38402	 0.22117	 0.11978
 1.000	 1.000	 0.500	 0.500	 1.000	 0.5000	 0.45552	 0.23091	 0.12708	 0.06238	 0.40992	 0.15991	 0.08006
 1.000	 1.000	 0.500	 0.500	 1.500	 0.5000	 0.62033	 0.17566	 0.09032	 0.03485	 0.35581	 0.11136	 0.05336
 1.000	 1.000	 0.500	 0.500	 2.000	 0.5000	 0.75747	 0.14136	 0.06675	 0.02105	 0.25822	 0.07372	 0.03713

- phi-dependence for U = 2:
 2.000	 1.000	 0.500	 0.500	 0.000	 0.0000	 0.26746	 0.50000	 0.18749	 0.26560	 0.26560	 0.26560	 0.00000
 2.000	 1.000	 0.500	 0.500	 0.000	 0.1000	 0.26045	 0.50000	 0.18744	 0.26464	 0.26464	 0.26464	 0.03044
 2.000	 1.000	 0.500	 0.500	 0.000	 0.3000	 0.20657	 0.50000	 0.18704	 0.25801	 0.25801	 0.25801	 0.08978
 2.000	 1.000	 0.500	 0.500	 0.000	 0.4000	 0.16216	 0.50000	 0.18677	 0.25342	 0.25342	 0.25342	 0.11827
 2.000	 1.000	 0.500	 0.500	 0.000	 0.5000	 0.10832	 0.50000	 0.18638	 0.24854	 0.24854	 0.24854	 0.14553
 2.000	 1.000	 0.500	 0.500	 0.000	 0.6000	 0.04666	 0.50000	 0.18592	 0.24368	 0.24368	 0.24368	 0.17122
 2.000	 1.000	 0.500	 0.500	 0.000	 0.6200	 0.03355	 0.50000	 0.18580	 0.24270	 0.24270	 0.24270	 0.17612
 2.000	 1.000	 0.500	 0.500	 0.000	 0.6400	 0.02017	 0.50000	 0.18571	 0.24177	 0.24177	 0.24177	 0.18098
 2.000	 1.000	 0.500	 0.500	 0.000	 0.6600	 0.00658	 0.50000	 0.18561	 0.24083	 0.24083	 0.24083	 0.18574
 2.000	 1.000	 0.500	 0.500	 0.000	 0.6650	 0.00315	 0.50000	 0.18558	 0.24059	 0.24059	 0.24059	 0.18691
 2.000	 1.000	 0.500	 0.500	 0.000	 0.6670	 0.00177	 0.50000	 0.18557	 0.24050	 0.24050	 0.24050	 0.18739
 2.000	 1.000	 0.500	 0.500	 0.000	 0.6680	 0.00108	 0.50000	 0.18556	 0.24045	 0.24045	 0.24045	 0.18762
 2.000	 1.000	 0.500	 0.500	 0.000	 0.6690	 0.00039	 0.50000	 0.18556	 0.24041	 0.24041	 0.24041	 0.18786
 2.000	 1.000	 0.500	 0.500	 0.000	 0.6695	 0.00004	 0.50000	 0.18556	 0.24040	 0.24040	 0.24040	 0.18798
```
Results may vary slightly, depending on the parameters set in *squad.in*.  

#### References:
1. [M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Phys. Rev. B 93, 024523 (2016).](http://dx.doi.org/10.1103/PhysRevB.93.024523)  
2. [M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Sci. Rep. 5, 8821 (2015).](http://dx.doi.org/10.1038/srep08821)  

