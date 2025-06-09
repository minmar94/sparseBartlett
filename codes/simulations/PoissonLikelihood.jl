# This script runs simulations for different graphical models using the SparseBartlettPackage in Julia.
# It generates synthetic count data (Poisson) under various settings (banded, random, sparse, etc.),
# introduces missingness and zero-inflation as specified, and runs MCMC for Bayesian inference using a GLM approach.
# The results are saved for further analysis in R.
#
# Key sections:
# - Parameter setup and environment activation
# - Data generation for different model structures
# - Handling of missing data and zero-inflation
# - MCMC sampling using the SparseBartlettPackage (GLM version)
# - Exporting results to R for post-processing

#### #### #### #### #### #### 
#### !Parameters
#### #### #### #### #### #### 

# Flag to indicate if this is the first run (for environment setup)
is_first_time::Bool = false # it is used to instantiate the environment
# use it only the first time youy ran any simulations

#### #### #### #### #### #### 
#### !simulations
#### #### #### #### #### #### 
# Activate the simulation environment and handle package setup if first run
using Pkg
Pkg.activate("./codes/simulations")

if is_first_time
    Pkg.rm("SparseBartlettPackage")
    Pkg.develop(url="SparseBartlettPackage")
    Pkg.instantiate()
end

# Load required Julia and R packages
using MKL
using RCall
using Distributions
using PDMats
using LinearAlgebra
using ToggleableAsserts
using Enzyme
using Random
using SparseBartlettPackage
## sim

# Load R BDgraph package if first run
if is_first_time
    RCall.@rlibrary("BDgraph")
end

# name of the output files
name_gen = "Simulation_"

# Main simulation loop over seeds and model sizes
n = 100
for i_seed in 1:20
    
    index_seed::Int64 = i_seed

    for i_k in 1:2
        
        k::Int64 = [10,25][i_k]
        
        # Loop over different data/model types. Data 4 in now shown in the paper and hence the results are not produced here
        for index_data in [1,2,3,5]

            # Set number of zero-inflation and decay scenarios based on model type
            if index_data == 5
                end_zero=4
                #end_pi=1
                end_decay=1
            else
                end_zero=1
                #end_pi=1
                end_decay=1
            end
            
            # Loop over missingness scenarios
            for index_miss in 1:2

                # Loop over zero-inflation scenarios
                for index_zero in 1:end_zero

                    # Loop over decay scenarios (for random graph types)
                    for index_decay in 1:end_decay

                        # Set likelihood and simulation parameters for the current scenario
                        index_like = 2
                        par_pi::Float64 = 0.5
                        par_like = [Normal(0.1, 0.1), Poisson(0.1), Bernoulli(0.1)][index_like]
                        par_miss = [0.0, 0.1][index_miss]
                        par_zero = [0, 0.25, 0.5, 0.75][index_zero]

                        # Generate a descriptive name for the simulation output
                        name = name_gen * "Model_" * string(["Banded1", "Banded2", "Banded3", "Random", "FullRandom"][index_data]) * "_Like" * string(["Normal", "Poisson", "Bernoulli"][index_like]) * "_pi=" * string(par_pi) * "_percmiss=" * string(par_miss)
                        name = name * "_n=" * string(n) * "_k=" * string(k)
                        if index_data == 4
                            name = name * "perczero=" * string(par_zero)
                            name = name * "decay=" * string([1.0, 3.0, 6.0][index_corr])
                        end
                        if index_data == 5
                            name = name * "perczero=" * string(par_zero)
                        end
                        ## sim
                        seed = 1
                        Random.seed!(seed + index_seed)
                        name = name * "seed=" * string(seed + index_seed)
                        #### sim
                        println(name)
                        println(typeof(par_like))
                        println([n, k])
                        kappa = 0.1
                        lambda = zeros(Float64, k, k)
                        z_mat = zeros(Int64, k, k)

                        y = zeros(Float64, k, n)::Matrix{Float64}
                        w = zeros(Float64, k, n)::Matrix{Float64}

                        beta_vec = rand(Normal(0.0, 2.0^0.5), k)
                        mu_mat = rand(Normal(0.0, 2.0^0.5), k, n)
                        cov = ones(Float64, k, 1, n)
                        
                        if index_data <= 3
                            # Banded models: set up banded precision and adjacency matrices
                            kappa = 1.0
                            for i = 1:k
                                lambda[i, i] = kappa
                                z_mat[i, i] = 1
                                for j = (i+1):min(k, i + index_data)
                                    lambda[i, j] = (-kappa * 1.0 + kappa * 0.001) / (index_data * 2.0)
                                    z_mat[j, i] = 1
                                    lambda[j, i] = (-kappa * 1.0 + kappa * 0.001) / (index_data * 2.0)
                                end
                            end
                            for i = 1:(k-1)
                                for j = (i+1):k
                                    z_mat[i, j] = 1
                                end
                            end
                        elseif index_data == 4
                            # Random graph: generate random covariance and precision matrices with spatial decay
                            xxx = rand(Uniform(0.0, 1.0), k)
                            yyy = rand(Uniform(0.0, 1.0), k)
                            decay = [1.0, 3.0, 6.0][index_corr]
                            varvv = 1.0
                            cov_mat = zeros(Float64, k, k)
                            for ik1 = 1:k
                                for ik2 = 1:k
                                    cov_mat[ik1, ik2] = varvv * exp(-decay * abs(sqrt((xxx[ik1] - xxx[ik2])^2 + (yyy[ik1] - yyy[ik2])^2)))
                                end
                            end
                            chol = cholesky(Symmetric(inv(cov_mat))).L
                            nm = Int64(((k^2 - k) / 2.0))
                            n_sparse::Int64 = Int64(trunc(nm * par_zero))
                            index_sparse::Vector{Int64} = sample(1:nm, n_sparse)
                            index_sparse = index_sparse[sortperm(index_sparse)]
                            index_not_sparse::Vector{Int64} = (1:nm)[(!in).((1:nm), Ref(index_sparse))]
                            row_col_indm::Matrix{Int64} = zeros(Int64, nm, 3)
                            h = 1
                            for j = 1:(k-1)
                                for i = (j+1):k
                                    row_col_indm[h, 1] = i
                                    row_col_indm[h, 2] = j
                                    row_col_indm[h, 3] = h
                                    h = h + 1
                                end
                            end
                            B = zeros(Float64, k, k)::Matrix{Float64}
                            vec_m = rand(Normal(0.0, 1.0), nm)::Vector{Float64}
                            vec_c = zeros(Float64, k)::Vector{Float64}
                            for ik = 1:k
                                vec_c[ik] = chol[ik, ik]
                            end
                            hhh = 0
                            for ik = 1:(k-1)
                                for ij = (ik+1):k
                                    hhh = hhh + 1
                                    vec_m[hhh] = chol[ij, ik]
                                end
                            end
                            nu = k + 1
                            h = 1
                            for j = 1:(k-1)
                                for i = (j+1):k
                                    B[i, j] = vec_m[h]
                                    h = h + 1
                                end
                                B[j, j] = vec_c[j]
                            end
                            B[k, k] = vec_c[k]
                            Q = B::Matrix{Float64}
                            lambda_withoutzero = Q * transpose(Q)
                            lambda_withoutzero
                            # set zeros in the precision matrix for sparsity
                            for h = index_sparse
                                i = row_col_indm[h, 1]
                                j = row_col_indm[h, 2]
                                if j == 1
                                    Q[i, j] = 0.0
                                else
                                    Q[i, j] = -sum(Q[i, 1:(j-1)] .* Q[j, 1:(j-1)]) / Q[j, j]
                                end
                            end
                            lambda .= Q * transpose(Q)
                            z_mat[:, :] .= 1
                            for i = 1:size(index_sparse, 1)
                                z_mat[row_col_indm[index_sparse[i], 1], row_col_indm[index_sparse[i], 2]] = 0
                            end
                            inv(lambda)
                        elseif index_data == 5
                            # Full random sparse graph: generate random sparse precision matrix using R's BDgraph
                            nm = Int64(((k^2 - k) / 2.0))
                            n_sparse = Int64(trunc(nm * par_zero))::Int64
                            index_sparse = sample(1:nm, n_sparse)::Vector{Int64}
                            index_sparse = index_sparse[sortperm(index_sparse)]
                            index_not_sparse = (1:nm)[(!in).((1:nm), Ref(index_sparse))]::Vector{Int64}
                            row_col_indm = zeros(Int64, nm, 3)::Matrix{Int64}
                            h = 1
                            for j = 1:(k-1)
                                for i = (j+1):k
                                    row_col_indm[h, 1] = i
                                    row_col_indm[h, 2] = j
                                    row_col_indm[h, 3] = h
                                    h = h + 1
                                end
                            end
                            z_mat[:, :] .= 1
                            for i = 1:size(index_sparse, 1)
                                z_mat[row_col_indm[index_sparse[i], 1], row_col_indm[index_sparse[i], 2]] = 0
                            end
                            @rput z_mat
                            @rput k
                            R"""
                                library(BDgraph)
                                mat = rgwish(1,t(z_mat),3,diag(1,k)) 
                            """
                            @rget mat
                            lambda .= mat
                            #inv(lambda)
                        else
                            error("T3")
                        end

                        # Standardize the precision matrix
                        app_sigma = inv(lambda)
                        app_sd = zeros(Float64, k, k)
                        for uuu = 1:k
                            app_sd[uuu, uuu] = app_sigma[uuu, uuu]^0.5
                        end
                        app_sd2 = (1.0 / (1.0^0.5)) .* Matrix(I(k))
                        lambda = app_sd2 * app_sd * lambda * app_sd * app_sd2

                        # Compute covariance matrix from precision matrix
                        sigma::Symmetric{Float64,Matrix{Float64}} = inv(Symmetric(lambda))
                        #beta_vec .= log.(rand(Uniform(2.0,7.0),k))
                        beta_vec = [log(5.0)]

                        # Generate mean structure and sample data from the Poisson GLM
                        for l = 1:n
                            mu_mat[:, l] = cov[:, :, l] * beta_vec
                            w[:, l] = rand(MvNormal(mu_mat[:, l], sigma))
                            for j = 1:k
                                y[j, l] = 1.0 * rand(Poisson(exp(w[j, l])))
                            end
                        end
                        y

                        # Handle missing data: randomly mask entries and impute based on likelihood
                        index_miss_mat = zeros(Int64, 0, 0)
                        n_miss = 0::Int64
                        index_tot = zeros(Int64, k * n, 2)
                        h = 0
                        true_y = zeros(Float64, 0)
                        for ii = 1:n
                            for ik = 1:k
                                h = h + 1
                                index_tot[h, 1] = ik
                                index_tot[h, 2] = ii
                            end
                        end
                        ydata = deepcopy(y)
                        if par_miss > 0.0
                            n_miss = Int64(trunc(ceil(n * k * par_miss)))
                            true_y = zeros(Float64, n_miss)
                            www = sample(1:(n*k), n_miss, replace=false)
                            sort!(www)
                            index_miss_mat = index_tot[www, :]
                            if typeof(par_like) == Normal{Float64}
                                for imiss = 1:n_miss
                                    true_y[imiss] = y[index_miss_mat[imiss, 1], index_miss_mat[imiss, 2]]
                                    ydata[index_miss_mat[imiss, 1], index_miss_mat[imiss, 2]] = rand(Normal(0.0, 1.0))
                                end
                            elseif typeof(par_like) == Poisson{Float64}
                                for imiss = 1:n_miss
                                    true_y[imiss] = y[index_miss_mat[imiss, 1], index_miss_mat[imiss, 2]]
                                    ydata[index_miss_mat[imiss, 1], index_miss_mat[imiss, 2]] = rand(Poisson(1.0))
                                end
                            elseif typeof(par_like) == Bernoulli{Float64}
                                for imiss = 1:n_miss
                                    true_y[imiss] = y[index_miss_mat[imiss, 1], index_miss_mat[imiss, 2]]
                                    ydata[index_miss_mat[imiss, 1], index_miss_mat[imiss, 2]] = rand(Bernoulli(0.5))
                                end
                            else
                                error("Like")
                            end
                        end

                        ##### ##### ##### ##### ##### ##### 
                        ##### MCMC 
                        ##### ##### ##### ##### ##### ##### 
                        z_mat_init = ones(Float64, k, k)

                        molt = 1
                        Random.seed!(seed + index_seed)
                        nm = Int64((k * k - k) / 2)
                        mcmc_par = (iter=10000 * molt, burnin=8000 * molt, thin=1 * molt)
                        mcmc_par = (iter=3, burnin=1, thin=1)
                        
                        out = mcmc_glm(
                            # Observed data matrix (with missing values imputed)
                            ydata=ydata,
                            # MCMC iteration parameters: (iter, burnin, thin)
                            iterations=mcmc_par,
                            # Covariate array for the mean structure
                            covariates=cov,
                            # Prior for the sparsity pattern of the precision matrix
                            model_sparse=GenSparse(Int64.(nm:(-1):1) ./ sum(Int64.(nm:(-1):1))),
                            # Prior for the scale/variance parameters
                            model_psi=GeneralS(),
                            # HMC sampler settings for the precision matrix
                            type_mcmc_lambda=HMC_DualAv_NoUTurn(; epsilon=10.000, L=1, M=1, restart=false, delta=0.5, iter_adapt=10000000),
                            # Prior for the precision matrix (Wishart distribution)
                            prior_lambda=Wishart(1.0 * (k + 3 - 1), 1.0 .* Matrix(I(k))),
                            # Prior for the mean parameter
                            prior_mu=Normal(0.0, 1000.0^0.5), 
                            # Initial values for scale parameters
                            c_init=ones(Float64, k),
                            # Initial values for mean vector
                            m_init=zeros(Float64, nm),
                            # Initial values for regression coefficients
                            beta_init=zeros(Float64, size(cov, 2)),
                            # Prior probability for edge inclusion in the graph
                            prob_zeta=ones(Float64, k, k) * par_pi,
                            # Initial adjacency matrix for the graph
                            z_init=Float64.(z_mat_init),
                            # Iteration to start HMC updates
                            iter_start_hmc=1,
                            # Iteration to start updating sparsity pattern
                            iter_start_zeta=500,
                            # Likelihood function (Poisson, Normal, Bernoulli)
                            like=par_like,
                            # Initial latent variable matrix for GLM
                            w_init=zeros(Float64, k, n),
                            # Indices of missing data in the observed matrix
                            missing_index=index_miss_mat
                        )
                        
                        # Extract MCMC output
                        m_out = out.m_out
                        c_out = out.c_out
                        qmat_sparse_out = out.q_mat_sparse_out
                        lambda_out = out.lambda_out
                        zmat_out = out.zmat_out
                        beta_out = out.beta_out
                        missing_out = out.missing_out
                        w_out = out.w_out
                        
                        # Export results and relevant variables to R for saving and further analysis
                        @rput true_y
                        @rput par_miss
                        @rput missing_out
                        @rput index_miss_mat
                        @rput beta_out
                        
                        @rput lambda_out
                        @rput name
                        @rput lambda
                        @rput k
                        @rput zmat_out
                        @rput z_mat
                        @rput beta_vec
                        @rput y
                        @rput w
                        @rput beta_vec
                        @rput w_out

                        dir_out = "./codes/simulations/out/"
                        @rput dir_out
                        R"""
                            save.image(paste(dir_out, name, "mcmc_out_lambda.Rdata", sep=""))
                        """
                        # End of script
                    end
                end
            end
        end
    end
end




