#### #### #### #### #### #### 
#### !Parameters
#### #### #### #### #### #### 

is_first_time::Bool = false # it is used to instantiate the environment
                            # use it only the first time youy ran any simulations

#### #### #### #### #### #### 
#### !simulations
#### #### #### #### #### #### 
using Pkg
Pkg.activate("./codes/simulations")


if is_first_time
    Pkg.rm("SparseBartlettPackage")
    Pkg.develop(url="SparseBartlettPackage")
    Pkg.instantiate()
end

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

if is_first_time

    RCall.@rlibrary("BDgraph")

end

n = 100
for i_seed in 1:20
    
    index_seed::Int64 = i_seed

    for i_k in [10, 25]
        
        k::Int64 = [10,25][i_k]
        
        for index_data in [1,2,3,5]

            if index_data == 5
                end_zero=4
                #end_pi=1
                end_decay=1
            else
                end_zero=1
                #end_pi=1
                end_decay=1
            end
            
            for index_miss in 1:2

                for index_zero in 1:end_zero

                    for index_decay in 1:end_decay

                    end

                end

            end
        end

    end

end
index_like = 1
par_pi::Float64 = 0.5
par_like = [Normal(0.1, 0.1), Poisson(0.1), Bernoulli(0.1)][index_like]
par_miss = [0.0, 0.1][index_miss]
par_zero = [0, 0.25, 0.5, 0.75][index_zero]

name = name_gen * "Model_" * string(["Banded1", "Banded2", "Banded3", "Random", "FullRandom"][index_data]) * "_Like" * string(["Normal", "Poisson", "Bernoulli"][index_like]) * "_pi=" * string(par_pi) * "_percmiss=" * string(par_miss)
name = name * "_n=" * string(n) * "_k=" * string(k)

if index_data == 5
    name = name * "perczero=" * string(par_zero)
end
## sim
seed = 1
Random.seed!(seed + index_seed);
name = name * "seed=" * string(seed + index_seed)
#### sim
lambda = zeros(Float64, k, k)
z_mat = zeros(Int64, k, k)

y::Matrix{Float64} = zeros(Float64, k, n);
w::Matrix{Float64} = zeros(Float64, k, n);

beta_vec = rand(Normal(0.0, 2.0^0.5), k)
mu_mat = rand(Normal(0.0, 2.0^0.5), k, n)
cov = ones(Float64, k, 1, n)


if index_data <= 3

    global kappa = 1.0
    for i = 1:k
        global lambda[i, i] = kappa
        global z_mat[i, i] = 1
        for j = (i+1):min(k, i + index_data)
            global lambda[i, j] = (-kappa * 1.0 + kappa * 0.001) / (index_data * 2.0)
            global z_mat[j, i] = 1
            global lambda[j, i] = (-kappa * 1.0 + kappa * 0.001) / (index_data * 2.0)
        end
    end



    for i = 1:(k-1)
        for j = (i+1):k
            global z_mat[i, j] = 1
        end
    end

elseif index_data == 4

    xxx = rand(Uniform(0.0, 1.0), k)
    yyy = rand(Uniform(0.0, 1.0), k)
    #decay = 3.0
    #decay = 1.0
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
            global h = h + 1
        end
    end



    B::Matrix{Float64} = zeros(Float64, k, k)
    vec_m::Vector{Float64} = rand(Normal(0.0, 1.0), nm)
    vec_c::Vector{Float64} = zeros(Float64, k)
    for ik = 1:k
        vec_c[ik] = chol[ik, ik]
    end
    hhh = 0
    for ik = 1:(k-1)
        for ij = (ik+1):k
            global hhh = hhh + 1
            vec_m[hhh] = chol[ij, ik]
        end

    end






    nu = k + 1
    
    global h = 1
    for j = 1:(k-1)
        for i = (j+1):k
            B[i, j] = vec_m[h]
            global h = h + 1
        end
        B[j, j] = vec_c[j]
    end
    B[k, k] = vec_c[k]
    Q::Matrix{Float64} = B
    lambda_withoutzero = Q * transpose(Q)
    lambda_withoutzero




    # set zeros
    for h = index_sparse
        i = row_col_indm[h, 1]
        j = row_col_indm[h, 2]
        if j == 1
            Q[i, j] = 0.0
        else
            Q[i, j] = -sum(Q[i, 1:(j-1)] .* Q[j, 1:(j-1)]) / Q[j, j]
        end
    end

    global lambda .= Q * transpose(Q)
    global z_mat[:, :] .= 1
    for i = 1:size(index_sparse, 1)
        global z_mat[row_col_indm[index_sparse[i], 1], row_col_indm[index_sparse[i], 2]] = 0
    end
    inv(lambda)

elseif index_data == 5



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
            global h = h + 1
        end
    end


    global z_mat[:, :] .= 1
    for i = 1:size(index_sparse, 1)
        global z_mat[row_col_indm[index_sparse[i], 1], row_col_indm[index_sparse[i], 2]] = 0
    end
    @rput z_mat
    @rput k
    R"""
        library(BDgraph)

        mat = rgwish(1,t(z_mat),3,diag(1,k))
    """
    @rget mat
    global lambda .= mat

    #inv(lambda)
else
    error("")
end

app_sigma = inv(lambda);
app_sd = zeros(Float64, k, k)
for uuu = 1:k
    app_sd[uuu, uuu] = app_sigma[uuu, uuu]^0.5
end
app_sd2 = (1.0 / (1.0^0.5)) .* Matrix(I(k))
lambda = app_sd2 * app_sd * lambda * app_sd * app_sd2



sigma::Symmetric{Float64,Matrix{Float64}} = inv(Symmetric(lambda))

beta_vec = [0.0]
for l = 1:n
    global mu_mat[:, l] = cov[:, :, l] * beta_vec
    y[:, l] = rand(MvNormal(mu_mat[:, l], sigma))
end


index_miss_mat = zeros(Int64, 0, 0)
n_miss::Int64 = 0
index_tot = zeros(Int64, k * n, 2)
h = 0
true_y = zeros(Float64, 0)
for ii = 1:n
    for ik = 1:k
        global h = h + 1
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
Random.seed!(seed + index_seed);
nm = Int64((k * k - k) / 2)
mcmc_par = (iter=10000 * molt, burnin=8000 * molt, thin=1 * molt)
Random.seed!(seed + index_seed);


out = mcmc(
    ydata=ydata,
    iterations=mcmc_par,
    covariates=cov,
    model_sparse=GenSparse(Int64.(nm:(-1):1) ./ sum(Int64.(nm:(-1):1))),
    model_psi=GeneralS(),
    type_mcmc_lambda=HMC_DualAv_NoUTurn(; epsilon=10.000, L=1, M=1, restart=false, delta=0.5, iter_adapt=10000000),
    prior_lambda=Wishart(1.0 * (k + 3 - 1), 1.0 .* Matrix(I(k))),
    prior_mu=Normal(0.0, 1000.0^0.5), c_init=ones(Float64, k),
    m_init=zeros(Float64, nm),
    beta_init=zeros(Float64, size(cov, 2)),
    prob_zeta=ones(Float64, k, k) * par_pi,
    z_init=Float64.(z_mat_init),
    iter_start_hmc=1,
    iter_start_zeta=500,
    missing_index=index_miss_mat
)

m_out = out.m_out;
c_out = out.c_out;
qmat_sparse_out = out.q_mat_sparse_out;
lambda_out = out.lambda_out;
zmat_out = out.zmat_out;
beta_out = out.beta_out;
missing_out = out.missing_out;

@rput elapsed;
@rput true_y;
@rput par_miss;
@rput missing_out;
@rput index_miss_mat;
@rput beta_out;

@rput lambda_out;
@rput dir_plot;
@rput name;
@rput lambda;
@rput k;
@rput zmat_out;
@rput z_mat;
@rput beta_vec;
@rput y;
@rput beta_vec;

R"""
    save.image(paste(dir_plot, name, "mcmc_out_lambda.Rdata", sep=""))
"""
