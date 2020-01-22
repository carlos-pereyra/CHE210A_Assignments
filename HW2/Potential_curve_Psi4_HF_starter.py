#%%
#            Chem 210A  2019  Prof. McCurdy
#
"""  Template for calculating a potential curve for a diatomic or one degree of freedom 
      for a polyatomic molecule

        using psi4 as a python module
  
    Starting point for Chem 210A problems -- CWM January 7, 2019

    This example does Hartree-Fock on H2
    
    Minor changes allow MP2 and CISD calculations on this and other diatomics
    With suitable modifications in the geometry specification 
    this template can be used to calculate a bond 
    stretching curve for a polyatomic
      
"""
#%%
#                    IMPORT MODULES
#
import psi4  #  import all of psi4
#  regular expression library is used for replacing strings with numbers in geometry specification
import re as re  # regular expression library 
# numpy library is used for all numerical calculations in python that we program specifically here 
import numpy as np # numpy library for numerical calculations
import os     # operating system commands for manipulating files
# plotting
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from matplotlib import rc
#
#                 SET MEMORY
#
psi4.set_memory('2048 MB')     # minimum memory should be 1 GB, 2GB or greater is better
#
#
#        SET FILE NAMES FOR OUTPUT FILES 
#         change names between runs to save them
#
#  File for output from Psi4 is a  .txt file 
file_string_output = 'output_potential_curve.txt'  # output file
#  molden output must be a .molden file
molden_file_name = "diatomic_scf.molden"   # molden file 
#  File for plotting potential curve 
plot_file_name = 'potcurve.dat'      # plotting file
#
# remove these files, especially the molden file if it exists
# to prevent psi4 from appending to existing file
if os.path.exists(molden_file_name):
  os.remove(molden_file_name)
if os.path.exists(plot_file_name):
  os.remove(plot_file_name)
if os.path.exists(file_string_output):
  os.remove(file_string_output)
psi4.core.set_output_file(file_string_output,False) 
print("psi4 output is directed to ",file_string_output,"\n")
#%%
#                    SPECIFY GEOMETRY PYTHON OBJECT
#
# R_value  =  the internuclear distance in a0
# and must be replaced by a number using the regular expression functions
# in the for loop on distances.
#
# Specifying cartesian coordinates in this example
# Types of atoms are specified by chemical symbol
# First line is charge and spin multiplicity 
# 0 1 => neutral  singlet 
# 0 3 => neutral triplet and requires ROHF ('rohf' reference for 'scf' below)
# c1 = C_1 symmetry (no symmetry) OK for RHF (and necessary for unrestricted Hartree-Fock (UHF))
# d2h => is fastest for homonuclear diatomics.  Irreps are (in this order in some commands for open shell systems)
#   Ag   B1g   B2g   B3g    Au   B1u   B2u   B3u 
# c2v => is fastest for heteronuclear diatomics.  Irreps are (in this order in some commands for open shell systems)#
#   A1    A2    B1    B2 
#
molecular_geometry  =  """
      0 1 
      H  0.0 0.0 0.0 
      F  0.0 0.0 R_value 
    units bohr
    symmetry c2v """

#   print the geometry template
print(molecular_geometry)
#
#  psi4 options that do not change with R value
# 
#                               BASIS SET
#  cc-pVDZ is a small correlation consistent Dunning basis - good for first calculation.
#  (Sometimes won't be  good enough for Chem 210A final answers
#  see Psi4 manual for many other built-in basis sets
#
psi4.set_options({'basis': 'cc-pvDZ'}) 
#
#                        METHOD (OR REFERENCE METHOD for MP2, CISD or CCSD)
#   Method, and/or  reference for CI or MP2 calculations is set here
# "rhf" for restricted HF for closed shells, 
# "rohf" for restricted open-shell HF for open shells, like triplets, requires other options to be set 
psi4.set_options({"reference": "rhf"})  
#  scf_type is a keyword that controls the algorithm for handling the four-index electron repulsion integrals.  
# 'pk' is a safe out-of-core algorithm using exact electron repulsion integrals.  
# 'df' denotes a density-fitted algorithm designed for computations with thousands of basis functions. 
#  It is required for MP2 (Moller-Plesset 2nd order perturbation theory)
psi4.set_options({'scf_type': 'pk'}) 
#
#                          SCF INITIAL GUESS
#  Initial guess options for SCF, very important for UHF and sometimes ROHF or RHF
#  RHF is usually not difficult to converge.
#
psi4.set_options({'guess': 'read'})  #  can be used to follow a solution for difficult cases as R changes
#                                 #  if we run an initial point before the loop on R with another option
#                                 #  then set guess: 'read' before starting the loop which is used for subsequent points
#   guess: 'core' good for closed shell restricted Hartree-Fock 
#psi4.set_options({'guess': 'core'})    # diagonalize core (nuclear attraction only) Hamiltonian for guess
#
#  some other guess options that help converge rohf and some rhf  calculations 
#psi4.set_options({'guess': 'sad'})    # superposition of atomic densities (SAD) guess
#psi4.set_options({'guess': 'gwh'})    # generalized Wolfsberg-Helmholz Huckel-like guess
psi4.set_options({"MAXITER": 500})    # If near 500 iterations are necessary, calculation is on the edge of failure
#%%
#                SINGLE POINT CALCULATION 
#      uncomment to do single point calculations 
#
#R = 1.73
#molecular_geometry_R = re.sub("R_value",str(R),molecular_geometry)
#scf_energy, wfcn = psi4.energy('scf',molecule=psi4.geometry(molecular_geometry_R),return_wfn=True)
#print("\n Single point calculation at R = ",R,"   psi4 scf energy = ", scf_energy)
#print(" Writing molden file ",molden_file_name)
#psi4.molden(wfcn,molden_file_name)
#exit()
#%%
#                    POTENTIAL CURVE CALCULATION
#
#  parameters for the curve
#
N_Rvals = 61  # when debugging use a small number of R values
Rmin = 1.0    # use a larger or smaller Rmin for some molecules (using bohr here) 
Rmax = 10.0
dr = (Rmax-Rmin)/(N_Rvals-1)
#
#  arrays for the R values and energies
#
Rvals=[]
Calculated_Energies=[]
#
#%%
#
# LOOP over R values 
#
print("\n ** Beginning potential curve calculation ** \n")
for n in range(N_Rvals):
    R = Rmin + n*dr
#
# The following command substitutes a string made from the value of the
# variable R for the the string "R_value" in the python object molecular_geometry
# to make a new object molecular_geometry_R .  It can be applied several times
# in succession, renaming the result each time, to replace several 
# strings to create a new geometry object
#
    molecular_geometry_R = re.sub("R_value",str(R),molecular_geometry)
#
#    print(molecular_geometry_R)  # print current geometry -- uncomment for debugging
#
#  method is specified here 'scf' assumes a reference to have been set
#  rhf,  rohf or uhf for example
#  'mp2'  also requires a reference (like rhf) as does 'cisd'
#
#               SCF calculation at this R point 
    scf_energy=psi4.energy('scf',molecule=psi4.geometry(molecular_geometry_R) )
    print("R = ",R," scf energy =  ", scf_energy)
    Rvals.append(R)
    Calculated_Energies.append(scf_energy)
#
#  for loop ends here
##
#%%
#                   MAKE OUTPUT FILE FOR PLOTTING
#  at the end of the for loop both print and make a file
#  containing the scf energies
#
f = open(plot_file_name,'w')  # optionally rename this file for the case at hand
print("\n Making plot file named ",plot_file_name)
for n in range(N_Rvals):
    print(Rvals[n],"   ",Calculated_Energies[n],file=f)
#   print(Rvals[n],"   ",Calculated_Energies[n]) # echo plotting output to console
#
#%%
#
#  MAKE PLOTS FROM THE OUTPUT FILE OR PUT PLOTTING COMMANDS HERE
#  FOR A PROPERLY LABELED GRAPH 
#  plotting with xmgrace is recommended, and easier than plotting
#  several cases from this python script
#
figure(num=None, figsize=(8, 5), dpi=70, facecolor='w', edgecolor='k')
plt.figure(1)
plt.plot( Rvals, Calculated_Energies, marker='+', markeredgecolor='black', color='black', markevery=1, linestyle='', label="V(R)")
plt.title(r'Ground State Potential R[1,10]', fontsize=15)
#plt.ylim(bottom=0, top=15)
#plt.ylim(bottom=0)
plt.xlabel(r'R [Bohr]', fontsize=12)
plt.ylabel(r'V(R) [Hart]', fontsize=12)
plt.legend(loc='best')
plt.show()
#
# problem 4. e.)
Ediff = min(Calculated_Energies)-Calculated_Energies[-1]
print("dE=E(min)-E(inf)={} [hartree]".format(Ediff))
print("dE=E(min)-E(inf)={} [kJ/mol]".format(Ediff*2625.5))
#
exit()

