# R2CageExperiment

Stochastic individual-based simulation and Approximate Bayesian Computation (ABC) fitting framework for a doublesex (*dsx*) gene-drive cage experiment in *Anopheles* mosquitoes that includes the R2 resistance allele.

The model tracks six genotypes (wild-type *W*, drive *D* / qfs5, and R2 resistance *R*) in both sexes and is used to match experimental pupal phenotype data (sex ratio and fluorescent marker frequencies).

## Repository Contents

| File / Directory | Description |
|------------------|-------------|
| `StochModel.cpp` + `head.h` | C++ stochastic generational model (source) |
| `CageSimulation` | Pre-compiled binary of the model (Linux/macOS-style executable) |
| `SimulationProcessingFunctions.wl` | Wolfram Language functions for phenotype sampling, distance metrics, ABC particle generation, and plotting |
| `RunModel.nb` | Main Mathematica notebook that sets parameters, runs ABC sampling, and generates posterior predictive plots |
| `Pupadata` | Experimental pupal count data (Wolfram expression) for 3 cages × 8 generations × 7 phenotypes |
| `ExPlots` | Pre-generated experimental data plots (imported by the processing functions) |
| `LICENSE` | MIT License |

## Model Overview

The simulation is a discrete-generation stochastic model of a closed cage population:

- Females and males of each genotype mate according to a mating probability.
- Mated females produce egg clutches with genotype-dependent fecundity (`phi` × relative fitness `omega`).
- Inheritance is governed by a 6×6×6 probability matrix that includes:
  - Homing rate (`e`)
  - Rate of non-homologous end joining (NHEJ) that generates R2 resistance alleles (`nj`)
- A fixed number of eggs (`gen_size`) is sampled to form the next generation (with optional larval mortality applied downstream).
- Output is the number of individuals of each of the 12 sex × genotype classes each generation.

Default parameter meanings (from `RunModel.nb`):

1. `phi` – wild-type egg clutch size  
2.  `omegaD` – relative female fitness of drive heterozygotes  
3. `omegaR` – relative female fitness of R2 heterozygotes  
4. `play` – probability a mated female lays a clutch  
5. `e` – homing rate  
6. `nj` – rate of R-allele creation via NHEJ  
7. `mate` – probability a female mates  
8. `gen_size` – number of eggs reared each generation  
9. `num_generations`  
10. `num_reps`  
+ initial counts of the 12 sex × genotype classes

## Requirements

- **Mathematica** (tested with version 15; the notebook uses standard Wolfram Language functions)
- A C++ compiler (only needed if you want to recompile `StochModel.cpp`)
- The pre-built `CageSimulation` binary must be executable and in the working directory when the notebook runs

## Usage

1. Clone the repository and `cd` into it.
2. Make the binary executable (if necessary):
   ```bash
   chmod +x CageSimulation
3. Open RunModel.nb in Mathematica.
    Evaluate the cells. The notebook:
        Loads experimental data and processing functions
        Sets default parameters and initial conditions
        Performs rejection ABC sampling (prior sampling → simulation → distance to data)
        Retains the best-fitting particles
        Generates posterior predictive plots comparing simulated vs. experimental female proportion, RFP, and GFP frequencies

You can also call the processing functions directly from a Mathematica session after loading SimulationProcessingFunctions.wl.
Data

    Pupadata contains the raw experimental phenotype counts used for fitting.
    Phenotypes are aggregated into categories corresponding to marker combinations and sex (see phenotype[] in SimulationProcessingFunctions.wl).

Citation / Context

This code accompanies analysis of a doublesex gene-drive cage experiment that includes an R2 resistance allele.

MIT License — see the LICENSE file.
Recompiling the Model (optional)

To recompile the C++ model from source:
bash

g++ -O3 -std=c++17 -o CageSimulation StochModel.cpp

The resulting binary expects a single line of space-separated numeric parameters on standard input (in the order defined in the notebook) and writes generation-by-generation counts to standard output.
