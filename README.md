# IAV_intracellular_model
MATLAB scripts for the stochastic simulation of intracellular replication and population-level modeling of Influenza A Virus (IAV).

## Overview

This pipeline simulates the intracellular replication dynamics of IAV using a stochastic modeling framework based on the Gillespie algorithm with the tau-leaping approximation. It is organized into two components: single-cell intracellular simulations and population-level models that use single-cell outputs as input.

## Requirements

- MATLAB R2019b or later (uses `writematrix`, `readmatrix`, and `parpool`)
- No additional toolboxes are required beyond the MATLAB base

## Repository Structure

```
IAV_intracellular_model/
├── Intracellular_model/     # Single-cell stochastic simulation scripts
│   ├── Central_IAVGillespie_single.m    # Run a single simulation (Figure 2B, S1)
│   ├── Central_IAVGillespie_multiple.m  # Run batch simulations across MOIs (Figure 2A, 3A-B)
│   ├── model_IAV_GillespieTauLeaping.m  # Core simulation engine
│   ├── Propensities.m       # Reaction propensity calculations
│   ├── ReactMOI.m           # Reaction assembly for multi-variant infections
│   ├── ReactProts.m         # Protein synthesis reactions
│   ├── ReactProtsGen.m      # Protein synthesis with genotype tracking
│   ├── ReactVir.m           # Virion assembly reactions
│   ├── Read_MOI.m           # Post-processing of simulation output
│   ├── RepCycles.m          # Replication cycle tracking
│   ├── dtCal.m              # Adaptive timestep calculation
│   ├── model_gillespie.m    # Exact Gillespie algorithm
│   ├── model_tauL.m         # Tau-leaping algorithm
│   ├── nRNPapprox.m         # Quasi-steady-state approximation for nRNP steps
│   ├── randsample_WithoutReplacement.m  # Sampling utility
│   └── Stock_matrix_MDCK.xlsx  # Reaction network and parameter definitions
└── Population_model/        # Population-level model scripts
    ├── Population_mode_data.m          # Compare simulations to experimental data (Figure 4B-C)
    ├── Population_model_Ne.m           # Estimate effective population size Ne (Figure 4A)
    ├── Population_model_MOIvFitness.m  # Fitness × MOI parameter sweep (Figure 5B)
    └── Population_model_selection.m    # Selection dynamics over generations (Figure 5A)
```

## Quick Start

### 1. Run a single simulation

Open MATLAB, navigate to the `Intracellular_model/` directory, and run:

```matlab
Central_IAVGillespie_single
```

This runs one stochastic simulation with default parameters (MOI = 10, IniV = 1, 12-hour timecourse). Results are printed to the workspace and used to generate plots equivalent to Figure 2B and S1.

### 2. Run batch simulations

```matlab
Central_IAVGillespie_multiple
```

This sweeps MOI values from 1–20 with 1,000 replicates each, using a parallel pool if available. Results are written to a timestamped output directory (e.g., `Result_03-17-2026_14-30-00/`).

**Note:** Batch runs are computationally intensive. Expect several hours of runtime on a standard workstation.

### 3. Run population-level models

Population model scripts require pre-computed single-cell results. Either:
- Run `Central_IAVGillespie_multiple.m` to generate your own, or
- Download the pre-computed dataset from OSF: https://osf.io/f9e6v/

Then update the `current` variable at the top of each population script to point to your results folder:

```matlab
current = strcat(pwd, '/Result_YOUR-TIMESTAMP-HERE');
```

Then run any of the four population model scripts from within the `Population_model/` directory.

## Parameter Configuration (Excel File)

The reaction network, kinetic parameters, and initial conditions are defined in `Stock_matrix_MDCK.xlsx`. This file must contain six sheets:

| Sheet | Contents |
|-------|----------|
| 1 | **Stoichiometric matrix** — how each reaction changes species copy numbers (rows = species, columns = reactions) |
| 2 | **Reactants matrix** — which species are reactants in each reaction (used to compute propensities) |
| 3 | **Rate constants matrix** — kinetic parameter values per reaction |
| 4 | **Parameter list** — maps parameter symbols to numerical values |
| 5 | **QSSA-eliminated reactions** — parameters for reactions removed under the Quasi-Steady-State Approximation |
| 6 | **Virion composition** — stoichiometry of proteins in a mature IAV virion |

> **Note:** Degradation reactions are auto-generated internally, but degradation rate constants must be present in Sheet 4.

## Output File Format

`Central_IAVGillespie_multiple.m` writes CSV files to a timestamped folder for each MOI value. The filename prefix indicates MOI (e.g., `MOI_5_`):

| File | Contents |
|------|----------|
| `MOI_N_Genomes_*.csv` | Genome segment counts per replicate; rows = segments × replicates, columns = MOI variants |
| `MOI_N_Fit_*.csv` | Row 1: virion output per replicate; Row 2: simulation end times |
| `MOI_N_Virs_G_*.csv` | Genome composition of assembled virions; rows = segments × max virions, columns = replicates |
| `MOI_N_Virs_C_*.csv` | Virion composition (protein counts); same shape as Virs_G |
| `MOI_N_Virs_time_*.csv` | Time of virion assembly for each virion per replicate |
| `MOI_N_Tit_*.csv` | Virion count timecourse; rows = timepoints, columns = replicates |
| `MOI_N_Times_*.csv` | Time values corresponding to Tit rows |

## Population Model Scripts

| Script | Description | Figure |
|--------|-------------|--------|
| `Population_mode_data.m` | Loads single-cell outputs and compares a neutral null model to IAV simulations. Tracks 2,048 barcoded variants and estimates Ne. Experimental data from [Varble et al. 2014](https://doi.org/10.1016/j.chom.2014.09.009) loaded from `Data.xlsx`. | 4B–C |
| `Population_model_Ne.m` | Estimates effective population size (Ne) by comparing variant frequencies before and after a single round of infection, under null and IAV models. | 4A |
| `Population_model_MOIvFitness.m` | Sweeps fitness advantage (1–5×) and MOI (1–10) in a 10×10 grid. Introduces one high-fitness variant and tracks frequency over 28 generations (~14 days). Results shown as a heatmap. | 5B |
| `Population_model_selection.m` | Simulates four selection scenarios varying MOI and fitness advantage. Tracks beneficial mutant frequency across multiple viral generations. | 5A |

## Citation

Segredo-Otero E & Gresham D (2025). *Stochastic intracellular replication dynamics shape population-level evolution in Influenza A Virus.* (in preparation)

Pre-computed simulation data: https://osf.io/f9e6v/
