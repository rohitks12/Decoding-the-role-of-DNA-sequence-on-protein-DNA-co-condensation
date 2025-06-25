# ESPResSo Simulation Files

This repository contains ESPResSo input files, type-ID mapping files, initial configuration files, and main simulation scripts for modeling DNA–protein co-condensation.

---

## Type ID Files  
**Files:** `TypeID_homo.dat`, `TypeID_hetero_I.dat`, `TypeID_hetero_II.dat`, `TypeID_lambda.dat`

Each file contains 2 columns and 500 rows:
The first column lists the particle IDs (or indices) corresponding to the 500 monomers in the polymer.
The second column specifies the type IDs assigned to each monomer.

These type IDs are used in the main simulation scripts:
`main_script_homogeneous.tcl`
`main_script_heterogeneous_I.tcl`
`main_script_heterogeneous_II.tcl`
`main_script_lambda.tcl`

They determine the monomer–protein interaction strengths during the simulation.

---

## Initial Configuration Files  
**Files:** `config_Re_0.2.dat`, `config_Re_0.4.dat`, `config_Re_0.6.dat`, `config_Re_0.8.dat`

These files provide the initial spatial configuration of the polymer chain and are used as input in the simulation scripts.  
Each file contains 500 rows and 4 columns:
Column 1: Particle ID (0 to 499),
Columns 2–4: Initial coordinates (x, y, z) of each monomer in the three-dimensional simulation box

---

## Simulation Scripts

`main_script_homogeneous.tcl`  
  Loads `TypeID_homo.dat`, initializes a homogeneous polymer, and sets interaction parameters for monomer–protein interactions.

`main_script_heterogeneous_I.tcl` and `main_script_heterogeneous_II.tcl` 
  Load the corresponding heterogeneous type-ID files to simulate polymers with two types of monomers.

 `main_script_lambda.tcl`
  Uses `TypeID_lambda.dat` to simulate protein condensation with partial lambda DNA composition.

---

## Function File

This Tcl script defines utility procedures for managing simulation data in ESPResSo. It includes:

1. `save_sim`:  
   Saves the full simulation state, including global variables, Tcl variables, particle data, interaction parameters, and bond information to a specified file.

2. `save_sim1`:  
   Saves a minimal simulation state, including only global variables, Tcl variables, and particle data.

3. `readData`:  
   Reads an input file line by line, returning its contents as a list where each element corresponds to one line of the file.

---

## Notes

We have used ESPResSo-3.3.1 for performing these simulations.
It is nessary that the type ID file and initial configuration file used in a simulation match the setup expected by the script.
Parameters like protein concentration,interaction strength and runtime can be modified by editing the corresponding sections in the `.tcl` scripts.

