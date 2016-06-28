---
output: html_document
---
SQUAD
=====
**Description:**  
**SQUAD** (**su**perconducting **qua**ntum **d**ot) is a small python code that allows to calculate spectral and transport properties of a single-level quantum dot connected to two BCS superconducting (sc) leads, optinally including third, non-sc lead.

**Basic features:**  
- based on diagrammatic perturbation theory in Coulomb interaction U in Nambu formalism  
- calculates properties within the 2nd order perturbation theory ~U^2  
- weak-coupling expansion, works well for small U/Gamma  
- real-energy axis representation  
- direct access to spectral functions and Andreev bound states frequencies  
- reliable description of the zero-phase and zero-pi boundary  

**List of files:**  
- squadlib1.py - library of general functions and the HF solver  
- squadlib2.py- library of functions for calculating 2nd order PT  
- second_sc.py - main code to calculate 2nd order PT results for a system with two sc electrodes  
This code uses self-consistency in the static part of self-energies only.  
- second_fullsc.py - modified second_sc code that incoroprates full self-consistent calculation of self-energies, 
including the dynamical part. Unstable in some regions of parameter space.  
- run_2nd.py - script to run second_sc.py for a range of parameters  
- run_2nd_full.py - script to run second_fullsc.py for a range of parameters  
- ssn_second.py - main code to calculate 2nd order PT results in a system with one normal and two
sc electrodes  
- ssnlib.py - library of functions for calculating ssn results  
- params.py - script to read ssn.in file  
- ssn.in - parameter file for ssn_second.py, read by params.py  
- run_ssn.py - script to run ssn_second.py for a range of parameters  

**References:**  
- M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Phys. Rev. B 93, 024523 (2016).  
- M. Žonda, V. Pokorný, V. Janiš and T. Novotný, Sci. Rep. 5, 8821 (2015).  
