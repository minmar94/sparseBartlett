# A New Hierarchical Distribution on Arbitrary Sparse Precision Matrices

This repository contains data and code to reproduce the analysis in the paper _"A New Hierarchical Distribution on Arbitrary Sparse Precision Matrices"_, by G. Mastrantonio, P. Alaimo Di Loro, and M. Mingione.

The work introduces a general strategy for defining distributions over the space of sparse symmetric positive definite matrices. In particular, it uses the Cholesky factorization of the precision matrix, imposing sparsity through constraints on its elements while preserving independence and avoiding the need for numerical evaluation of normalization constants. The proposed prior is named **S-Bartlett**.

---

## Simulations

The `sparseBartlett/simulations` folder contains Julia scripts for running simulation studies using the `SparseBartlettPackage`. Each script generates synthetic data under different settings and runs MCMC for Bayesian inference.

### Prerequisites
- Julia (v1.9 or later recommended)
- All required Julia packages (see `Project.toml`)
- R (for scripts using `RCall` and `BDgraph`)
- The `SparseBartlettPackage` must be available in the workspace

### Running Simulations

To run the simulation with a normal likelihood:
```julia
include("simulations/NormalLikelihood.jl")
```

For the Poisson likelihood:
```julia
include("simulations/PoissonLikelihood.jl")
```

#### Notes
- If running a script for the first time, set `is_first_time = true` at the top of the file to install and instantiate dependencies. Set it to `false` for subsequent runs.
- Output files and results will be saved in the designated output directories (see the `dir_plot` or `dir_out` variables in each script).
- Some scripts require R and the `BDgraph` package. Make sure R is installed and accessible from Julia.

### Creating RData Files for Figures

The scripts `for_figures_normals.R` and `for_figures_poisson.R` (in the `simulations` folder) create the `.Rdata` files needed to produce the figures in the paper. Run these after completing the simulations.

> 💬 **MINGIOOOOO**: You might want to describe here how to generate the plots from the files in the `for_figures_` folders. I’ve left a couple of files for reference. Your figure data is already in the folder.

---

## Real Data Applications

The `applications` folder contains example scripts showing how to apply the `SparseBartlettPackage` to real datasets and statistical models.

### Prerequisites
- Julia (v1.9 or later recommended)
- All required Julia packages (see `Project.toml`)
- R (for some scripts using `RCall`)
- The `SparseBartlettPackage` must be available in the workspace

### Running the Models

To run the **Doubs** dataset, assuming independence between species:
```julia
include("applications/doubs.jl")
```


For the **genetic data** examples:
```julia
include("applications/geneExp50.jl")
include("applications/geneExp50_missing.jl")
include("applications/geneExp100.jl")
include("applications/geneExp100_missing.jl")
```
Choose the script based on the number of observations (50 or 100), and whether to include missing values for CRPS computation.

To run the **G-Wishart** model on the genetic data:
```julia
include("applications/gwishart_geneExp.R")
```

#### Notes
- As above, if running for the first time, set `is_first_time = true` at the top of the script.
- Output files and results will be saved in the specified output directory.
- Ensure R is installed and callable from Julia.
