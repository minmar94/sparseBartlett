name = "doubs_indipendent_species"

#### PARAMETERS ####
# Flag to indicate if this is the first run (for environment setup)
is_first_time::Bool = false # Set to true only the first time you run any simulations

#### SIMULATION ENVIRONMENT SETUP ####
using Pkg
Pkg.activate("./applications") # Activate the Julia environment for the applications folder
if is_first_time
    Pkg.rm("SparseBartlettPackage") # Remove the package if it exists (clean install)
    Pkg.develop(url="SparseBartlettPackage") # Develop the SparseBartlettPackage from local or remote
    Pkg.instantiate() # Install all dependencies
end

# Load required Julia packages
using MKL
using RCall
using Distributions
using PDMats
using LinearAlgebra
using ToggleableAsserts
using Enzyme
using Random
using SparseBartlettPackage

Random.seed!(1); # Set random seed for reproducibility

#### DATA LOADING AND PREPARATION ####
# Use R to load and preprocess the dataset
R"""
dataset = read.csv("./applications/data/doubs.csv", header = F)
summary(dataset)
site = dataset[1,]
species = dataset[,1]
dataset = dataset[-1,]
dataset = dataset[,-1]
dataset = as.matrix(dataset)
dataset <- t(dataset[,-8])
"""
@rget dataset; # Retrieve the processed dataset from R to Julia
dataset = 1.0 .* dataset # Ensure dataset is Float64
n::Int64 = size(dataset, 2) # Number of samples
k::Int64 = size(dataset, 1) # Number of species (variables)

#### MCMC SETUP ####
z_mat_init = ones(Float64, k, k) # Initial value for the sparsity pattern
molt = 1 # Multiplier for iterations (can be used to scale up/down)
nm = Int64((k * k - k) / 2) # Number of unique off-diagonal elements in a symmetric matrix

# Set MCMC parameters: number of iterations, burn-in, and thinning
mcmc_par = (iter=10000 * molt, burnin=8000 * molt, thin=1 * molt)

covariate = ones(Float64, k, 1, n) # Covariate matrix (all ones, no covariate effect)

#### RUN MCMC SAMPLER ####
out = mcmc_glm(
    ydata=dataset, # Observed data
    iterations=mcmc_par, # MCMC parameters
    covariates=covariate, # Covariate matrix
    model_sparse=GenSparse(Int64.(nm:(-1):1) ./ sum(Int64.(nm:(-1):1))), # Prior for sparsity pattern
    model_psi=GeneralS(), # Prior for psi parameter
    type_mcmc_lambda=HMC_DualAv_NoUTurn(; epsilon=10.000, L=10, M=10, restart=false, delta=0.5, iter_adapt=10000000), # HMC sampler settings
    prior_lambda=Wishart(1.0 * (k + 3 - 1), 1.0 .* Matrix(I(k))), # Wishart prior for precision matrix
    
    prior_mu=Normal(0.0, 1000.0^0.5), # Prior for mean
    c_init=ones(Float64, k), # Initial value for c parameter
    m_init=zeros(Float64, nm), # Initial value for m parameter
    beta_init=zeros(Float64, size(covariate, 2)), # Initial value for regression coefficients
    prob_zeta=ones(Float64, k, k) * 0.5, # Initial probability for zeta (sparsity)
    
    z_init=Float64.(z_mat_init), # Initial value for z (sparsity pattern)
    w_init=zeros(Float64, k, n), # Initial value for latent variables
    like=Poisson(1.0), # Poisson likelihood
    iter_start_hmc=2, # When to start HMC updates
    iter_start_zeta=5, # When to start zeta updates
    
)

#### EXTRACT AND SAVE RESULTS ####
# Extract MCMC outputs
m_out = out.m_out;
c_out = out.c_out;
qmat_sparse_out = out.q_mat_sparse_out;
lambda_out = out.lambda_out;
zmat_out = out.zmat_out;
beta_out = out.beta_out;
y = dataset;

# Send results to R for saving
@rput y;
@rput beta_out;
@rput lambda_out;
@rput name;
@rput k;
@rput zmat_out;

# Save the R workspace with all results for later analysis
R"""
    save.image(paste("./applications/out/",name, "mcmc_out_lambda.Rdata",sep=""))
"""
