---
output: pdf_document
---
SQUAD
=====
**Description:**  
**SQUAD** (**su**perconducting **qua**ntum **d**ot) is a set of small python codes that allow to calculate spectral 
and transport properties of a single-level quantum dot connected to two BCS superconducting (sc) leads, 
optionally including third, non-sc lead (codes for calculations with normal electrode have *ssn* in their names).
  
SQUAD requires SciPy libraries (www.scipy.org) and it was developed and tested using python 2.7.5 and SciPy 0.12.  
  
These codes are subject to big changes and require some optimization and cleanup. Please use with caution. 
In case of interest please contact me at pokornyv@fzu.cz.

[//]: # (![Setup with two superconducting and one optional normal electrode connected to correlated quantum dot](dot_lead_ssn2.jpg))

**Basic features:**  
- based on diagrammatic perturbation theory in Coulomb interaction U in Nambu formalism  
- calculates properties within the 2nd order perturbation theory ~U^2  
- weak-coupling expansion, works well for small U/Gamma  
- real-energy axis representation  
- direct access to spectral functions and Andreev bound states frequencies  
- reliable description of the zero-phase and zero-pi boundary  

**Solution for infinite Delta:**  
This code also allows to calculate an approximation of the system with infinitely large gap. Setting Delta large 
(e.g. 1e5 and larger) usually gives a solution numerically correct to 4-5 decimal places. Be sure you set Delta larger
thatn the maximum of the energy axis to get rid of the band contributions.  

**List of files:**  
- *squadlib1.py* - library of general functions and the HF solver  
- *squadlib2.py* - library of functions for calculating 2nd order PT  
- *second_sc.py* - main code to calculate 2nd order PT results for a system with two sc electrodes  
This code uses self-consistency in the static part of self-energies only.  
- *second_fullsc.py* - modified second_sc code that incoroprates full self-consistent calculation of self-energies, 
including the dynamical part. Unstable in some regions of parameter space.  
- *run_2nd.py* - script to run *second_sc.py* for a range of parameters  
- *run_2nd_full.py* - script to run *second_fullsc.py* for a range of parameters  
- *ssn_second.py* - main code to calculate 2nd order PT results in a system with one normal and two
sc electrodes  
- *ssnlib.py* - library of functions for calculating ssn results  
- *params.py* - script to read ssn.in file  
- *ssn.in* - parameter file for *ssn_second.py*, read by *params.py*, described in *ssn_infile.md*
- *run_ssn.py* - script to run *ssn_second.py* for a range of parameters  

**References:**  
- M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Phys. Rev. B 93, 024523 (2016).  
- M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Sci. Rep. 5, 8821 (2015).  
