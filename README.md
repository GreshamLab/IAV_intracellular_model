# IAV_intracellular_model
Matlab scripts for the stochastic simulation of intracellular replication and population level modeling.

================================================================================================================
IAV Stochastic Simulation Pipeline – Overview
================================================================================================================

This pipeline simulates the intracellular replication dynamics of Influenza A Virus (IAV) using a stochastic 
modeling framework based on the Gillespie algorithm with the tau-leaping approximation.

The pipeline is organized into two main MATLAB scripts:

- `Central_IAVGillespie_single.m`: Runs a single stochastic simulation for a specified MOI and initial condition.
- `Central_IAVGillespie_multiple.m`: Performs multiple replicate simulations for a range of MOI values or parameters,
  storing the aggregated outputs for downstream analysis or plotting.

These scripts call a core simulation module (`mode_IAV_GillespieTauLeaping.m`) that implements the tau-leaping 
method to efficiently simulate reaction dynamics in a stochastic setting.

------------------------------------------------------------------------------------------------
The model structure (reaction network, parameters, and initial conditions) is provided externally
via an Excel file (`Param_mat.xlsx`) with the following required sheets:

1. Stoichiometric matrix: Defines how species counts change per reaction.
2. Reactants matrix: Indicates which species are reactants in each reaction.
3. Rate constants matrix: Contains numeric rate parameters associated with each reaction.
4. Parameter list: Maps parameter names to values for interpretability and use in the simulation.
5. QSSA-eliminated reactions: Lists reactions removed under the Quasi-Steady-State Approximation (QSSA).
6. Virion composition: Specifies the stoichiometry of viral proteins in mature virions.

Note:
- Species-specific degradation reactions are automatically generated internally,
  but degradation rates must be defined in the parameter list.

------------------------------------------------------------------------------------------------
On execution, the scripts prompt the user (or accept input variables) for two key parameters:

- MOI: Multiplicity of Infection – the number of distinct viral variants infecting the cell.
- IniV: Initial number of virions per variant.

The total number of infecting particles is given by MOI × IniV.

------------------------------------------------------------------------------------------------
This framework enables the simulation of viral diversity, replication kinetics, and stochastic effects
within infected cells under varying initial conditions and parameter settings.

For example applications, see:
- `simulate_multiple.m` → used to generate data for Figure 2B and Figure S1 in Segredo-Otero & Gresham (2025)

