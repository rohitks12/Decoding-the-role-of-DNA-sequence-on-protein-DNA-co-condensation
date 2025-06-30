# ESPResSo Simulation Files

This repository contains ESPResSo input files, type ID files, initial configuration files, and main simulation scripts for modeling DNA–protein co-condensation.

---

## Type ID Files  
**Files:** `TypeID_homo.dat`, `TypeID_hetero_I.dat`, `TypeID_hetero_II.dat`, `TypeID_lambda.dat`

Each file contains 2 columns and 500 rows.
The first column lists the particle IDs (or indices) corresponding to the 500 monomers in the polymer.
The second column specifies the type IDs assigned to each monomer.

These type IDs are used in the main simulation scripts:
`main_script_homogeneous.tcl`,
`main_script_heterogeneous_I.tcl`,
`main_script_heterogeneous_II.tcl`, and
`main_script_lambda.tcl` 

They determine the monomer–protein interaction strengths for each monomer during the simulation, therefore define the sequence heterogeneity of DNA.


---

## Initial Configuration Files  
**Files:** `config_Re_0.2.dat`, `config_Re_0.4.dat`, `config_Re_0.6.dat`, `config_Re_0.8.dat`

These are the initial polymer configuration files used as inputs in the main simulation script. Each file contains 500 rows and 4 columns. The first column lists the particle IDs (or indices) corresponding to each monomer in the polymer. The last 3 columns specify the initial coordinates of these monomers in the three dimensional simulation box. 

---

## Simulation Scripts

`main_script_homogeneous.tcl`  
The file loads `TypeID_homo.dat` and initializes a homogeneous polymer, and sets interaction parameters for monomer–protein interactions along with other simulation parameters.

`main_script_heterogeneous_I.tcl` and `main_script_heterogeneous_II.tcl`  
The files load the corresponding heterogeneous type ID files to simulate polymers with two types of monomers.

`main_script_lambda.tcl`  
The file uses `TypeID_lambda.dat` to simulate protein condensation with partial lambda DNA composition.

---

## `functions.tcl`

This Tcl code defines procedures for managing simulation data in ESPResSo. It contains multiple functions.

1. `save_sim`: It saves the simulation state which includes global variables, tcl variables, particle data, interactions and bond informations to a specific file.
   
2. `save_sim1`: It saves only core simulation components such as variables, Tcl variables and particle data.

3. `readData`: It reads an input file line-by-line, returning its contents as a list where each element represents one line from the original file. 
---

## `force_calculation.py`

This python code computes the capillary forces exerted on "bare" DNA (DNA outside the condensate) by protein-DNA condensation. It uses the DBSCAN clustering algorithm to identify the DNA inside the condensate. Next, it computes the capillary force by calculating the bond stretch between the monomers present outside the condensate. 

## `probability_calculation.py`

This python code calculates the probability of each DNA monomer being part of the condensate. It first identifies monomers in the condensate for each output configuration, scoring them as 1 if they are inside the condensate and 0 if they are not. Next, it processes 3000 output files and averages these scores for each monomer to produce the probability profile. 

## `interfacial_affinity_calculation.py`

The code uses previously calculated probability profiles to identify interfacial monomers. The algorthm searches for increasing trends (left interface) and decreasing trends (right interface) in the probability profile. Once, the interfaces are identified, it computes the interfacial affinity at each interface and then averages these values to get the final average interfacial affinity.  

## `kymograph.py`

This python code generates a kymograph of DNA and protein distributions along the Z-axis over time. It extracts the Z-position of all the particle from output configurations. Next, it bins the data into equal spatial intervals, and normalizes by total particle count. Finally, it plots a heatmap with normalized end-to-end distance on horizontal axis and time on the vertical axis.  

## Notes

We have used ESPResSo-3.3.1 for performing these simulations.
It is necessary that the type ID file and initial configuration file used in a simulation match the setup expected by the script.
Parameters like protein concentration, interaction strength and runtime can be modified by editing the corresponding sections in the `.tcl` scripts.

