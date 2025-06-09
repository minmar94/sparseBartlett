

dir = "/home/gmastrantonio/Works/sparse/"
dir_plot = "/home/gmastrantonio/Works/sparse/out/"

name_gen = "Sim"
using Pkg
Pkg.activate(dir)


using MKL
using SparseIW
using RCall
using Distributions
using PDMats
using LinearAlgebra
using ToggleableAsserts
using Enzyme
using Random
## sim


n = parse(Int, ARGS[1])
k = [10,25,50][parse(Int, ARGS[2])]
index_pi = parse(Int, ARGS[3])
index_miss = parse(Int, ARGS[4]) 
index_data = parse(Int, ARGS[5])
index_zero = parse(Int, ARGS[6])
index_seed = parse(Int, ARGS[7])
index_corr = parse(Int, ARGS[8])

#n = 100
#k = 10
#index_pi = 4
index_like = 2
#index_miss = 2
#index_data = 3
#index_zero = 1

par_pi = [0.5, 0.1 , 0.9][index_pi]
#par_pi = [0.5][index_pi]
par_like = [Normal(0.1, 0.1), Poisson(0.1), Bernoulli(0.1)][index_like]
par_miss = [0.0,0.1][index_miss]
#par_miss = [0.0][index_miss]
par_zero = [0.1,0.25,0.5,0.75, 0.9][index_zero]

name = name_gen * "Model_" * string(["Banded1", "Banded2", "Banded3", "Random"][index_data]) * "_Like" * string(["Normal", "Poisson", "Bernoulli"][index_like]) * "_pi=" * string(par_pi) * "_percmiss=" * string(par_miss) 
name = name * "_n=" * string(n) * "_k=" * string(k)
if index_data == 4

    name = name * "perczero=" * string(par_zero)
    name = name * "decay=" * string([1.0,3.0,6.0][index_corr])
    

end


## sim
seed = 1
Random.seed!(seed+index_seed);  
name = name * "seed=" * string(seed+index_seed)
#### sim
println(name)
println(typeof(par_like))
println([n,k])
kappa = 0.1
lambda = zeros(Float64,k,k)
z_mat = zeros(Int64, k,k)

y::Matrix{Float64} = zeros(Float64, k, n);
w::Matrix{Float64} = zeros(Float64, k, n);

beta_vec = rand(Normal(0.0, 2.0^0.5),k)
mu_mat = rand(Normal(0.0, 2.0^0.5),k,n)
cov = ones(Float64,k,1,n)
#for iobs = 1:n
#    cov[:,:,iobs] .= 0.0
#    for ik =1:k
#        cov[ik,ik, iobs] = 1.0
#    end
#end
if index_data <= 3
    global kappa = 1.0
    for i = 1:k
        global lambda[i,i] = kappa
        global z_mat[i,i] = 1
        for j = (i+1):min(k,i+index_data)
            global lambda[i,j] = (-kappa *1.0 + kappa*0.001)/(index_data*2.0) 
            global z_mat[j,i] = 1
            global lambda[j,i] = (-kappa *1.0 + kappa*0.001)/(index_data*2.0)
        end
    end
    
    for i = 1:(k-1)
        for j = (i+1):k
            global z_mat[i,j] = 1
        end
    end
    


elseif index_data == 4
    xxx = rand(Uniform(0.0,1.0),k)
    yyy = rand(Uniform(0.0,1.0),k)
    ##decay = 3.0
 #decay = 1.0
 decay = [1.0,3.0,6.0][index_corr]
    varvv = 1.0
    cov_mat = zeros(Float64, k,k)
    for ik1 = 1:k

        for ik2 = 1:k
            cov_mat[ik1,ik2] = varvv* exp( - decay*abs(sqrt( (xxx[ik1] - xxx[ik2])^2 +  (yyy[ik1] - yyy[ik2])^2 )  )  )
        end
    end
    chol = cholesky(Symmetric(inv(cov_mat))).L
    
    #nu_psi = k+1
    #diag_psi = 2.0*(10.0*(nu_psi/2.0)-1.0)
    #app::Matrix{Float64} = zeros(Float64,k,k)
    #for i = 1:k
    #    app[i,i] = diag_psi
    #end
    #app2::Matrix{Float64} = zeros(Float64,k,k)
    #for i = 1:k
    #    app2[i,i] = 2.0
    #end
    ##psi::Matrix{Float64} = cholesky(Symmetric(rand(InverseWishart( nu_psi, app)))).L
    #psi = app2
    #psi = psi.*2.0
    #psi*psi'
    nm = Int64(((k^2-k)/2.0))

    n_sparse::Int64 = Int64(trunc(nm*par_zero))
    index_sparse::Vector{Int64} = sample(1:nm,n_sparse)
    index_sparse = index_sparse[sortperm(index_sparse)]

    index_not_sparse::Vector{Int64} = (1:nm)[(!in).((1:nm),Ref(index_sparse))]
    row_col_indm::Matrix{Int64} = zeros(Int64,nm,3)
    h = 1
    for j = 1:(k-1)
        for i = (j+1):k
            row_col_indm[h,1] = i
            row_col_indm[h,2] = j
            row_col_indm[h,3] = h
            global h = h+1
        end
    end 



    B::Matrix{Float64} = zeros(Float64,k,k)
    vec_m::Vector{Float64} = rand(Normal(0.0,1.0), nm)
    vec_c::Vector{Float64} = zeros(Float64,k)
    for ik = 1:k
        vec_c[ik] = chol[ik, ik]
    end
    hhh = 0
    for ik = 1:(k-1)
        for ij = (ik+1):k
            global hhh = hhh+1
            vec_m[hhh] = chol[ij, ik]
        end
        
    end

    


    #nu-1+1  = nu
    #nu-2+1 = nu-1


    nu  = k+1
    #for j = 1:k
    #    vec_c[j] = rand(Chisq(nu-(j-1-sum(row_col_indm[index_sparse ,1] .== j)) ))^0.5
    #    vec_c[j] = rand(Chisq(nu-(j-1) ))^0.5
    #end 
    global h = 1
    for j = 1:(k-1)
        for i = (j+1):k
            B[i,j] = vec_m[h]
            global h = h+1
        end
        B[j,j] = vec_c[j]
    end 
    B[k,k] = vec_c[k];
    Q::Matrix{Float64} = B
    lambda_withoutzero = Q*transpose(Q);
    lambda_withoutzero




    # set zeros
    for h = index_sparse
        i = row_col_indm[h,1]
        j = row_col_indm[h,2]
        if j == 1
            Q[i,j] = 0.0
        else
            Q[i,j] = - sum(Q[i,1:(j-1)].*Q[j,1:(j-1)])/Q[j,j]
        end
    end

    global lambda .= Q*transpose(Q);
    global z_mat[:,:] .= 1
    for i = 1:size(index_sparse, 1)
        global z_mat[row_col_indm[index_sparse[i], 1], row_col_indm[index_sparse[i], 2]] = 0
    end
    inv(lambda)

    #lambda .= lambda.*5.0
else 
    error("T3")
end

app_sigma = inv(lambda);
app_sd = zeros(Float64,k,k)
for uuu = 1:k
    app_sd[uuu,uuu] = app_sigma[uuu,uuu]^0.5
end
app_sd2 = (1.0/(1.0^0.5)).* Matrix(I(k))
lambda = app_sd2*app_sd*lambda*app_sd*app_sd2

#cholesky(lambda)
sigma::Symmetric{Float64,Matrix{Float64}} = inv(Symmetric(lambda))
#global beta_vec .= log.(rand(Uniform(2.0,7.0),k))
beta_vec = [log(5.0)]

for l = 1:n
    global mu_mat[:,l] = cov[:,:,l]*beta_vec
    w[:, l] = rand(MvNormal(mu_mat[:,l], sigma))
    for j = 1:k
        y[j, l] = 1.0*rand(Poisson(exp(w[j, l])))
    end
end
y

index_miss_mat = zeros(Int64,0,0)
n_miss::Int64 = 0
index_tot = zeros(Int64,k*n,2)
h = 0
true_y =  zeros(Float64,0)
for ii = 1:n
    for ik = 1:k
        global h = h+1
        index_tot[h,1] = ik
        index_tot[h,2] = ii
    end
end
ydata = deepcopy(y)
if par_miss > 0.0
    n_miss = Int64(trunc(ceil(n*k*par_miss)))

    true_y = zeros(Float64,n_miss)
    www = sample(1:(n*k), n_miss, replace=false)
    sort!(www)
    index_miss_mat = index_tot[www,:]
    if typeof(par_like) == Normal{Float64}
        for imiss = 1:n_miss
            true_y[imiss] = y[index_miss_mat[imiss,1], index_miss_mat[imiss,2]]
            ydata[index_miss_mat[imiss,1], index_miss_mat[imiss,2]] = rand(Normal(0.0,1.0))
        end
    elseif typeof(par_like) == Poisson{Float64}
        for imiss = 1:n_miss
            true_y[imiss] = y[index_miss_mat[imiss,1], index_miss_mat[imiss,2]]
            ydata[index_miss_mat[imiss,1], index_miss_mat[imiss,2]] = rand(Poisson(1.0))
        end
    elseif typeof(par_like) == Bernoulli{Float64}
        for imiss = 1:n_miss
            true_y[imiss] = y[index_miss_mat[imiss,1], index_miss_mat[imiss,2]]
            ydata[index_miss_mat[imiss,1], index_miss_mat[imiss,2]] = rand(Bernoulli(0.5))
        end
    else
        error("Like")
    end

end




##### ##### ##### ##### ##### ##### 
##### MCMC 
##### ##### ##### ##### ##### ##### 
z_mat_init = ones(Float64,k,k)

molt = 1
Random.seed!(seed+index_seed);  
nm = Int64((k*k- k)/2)
mcmc_par = (iter = 10000 * molt, burnin = 8000 * molt, thin = 1 * molt)
#mcmc_par = (iter = 100 , burnin = 80 , thin = 1  )
elapsed = @elapsed begin
out = mcmc_glm(
        ydata = ydata,
        iterations = mcmc_par,
        covariates = cov,
        #iterations = (iter = 100 , burnin = 20 , thin = 1 ),
        #model_sparse = NoSparse(),
        #model_sparse = GenSparse(reverse(pdf(Exponential(25.0), 0:nm)/sum(pdf(Exponential(25.0), 0:nm)))),
        #model_sparse = GenSparse(ones(Float64, nm+1)/(nm+1)),
        model_sparse = GenSparse(Int64.(nm:(-1):1)./sum(Int64.(nm:(-1):1))),
        model_psi = GeneralS(),
        #type_mcmc_lambda = BaseHMC(;epsilon = 1.0, L = 10, M = 10, restart = true),
        #type_mcmc_lambda = HMC_DualAv(;epsilon = 1.0, L = 10, M = 10, restart = false, delta = 0.5, lambda = 0.001, iter_adapt = 10000000   ),
        type_mcmc_lambda = HMC_DualAv_NoUTurn(;epsilon = 10.000, L = 10, M = 10, restart = false, delta = 0.5,  iter_adapt = 10000000   ),
        #type_mcmc = HMC(0.8),

        prior_lambda = Wishart(1.0*(k+3-1), 1.0.*Matrix(I(k))),
        prior_mu = Normal(0.0,1000.0^0.5),

        c_init = ones(Float64, k),
        m_init = zeros(Float64,nm),
        beta_init = zeros(Float64, size(cov,2)),
        prob_zeta = ones(Float64,k,k)*par_pi, 
        #m_init = vec_m,
        #c_init = vec_c,#[sqrt(rand(Chi((k+2)-i+1))) for i = 1:k ],
        z_init = Float64.(z_mat_init),
        iter_start_hmc = 2,
        iter_start_zeta = 500,
        like = par_like,
        w_init = zeros(Float64,k,n),
        missing_index = index_miss_mat
    );
end

m_out= out.m_out;
c_out = out.c_out;
qmat_sparse_out = out.q_mat_sparse_out;
lambda_out = out.lambda_out;
zmat_out = out.zmat_out;
beta_out = out.beta_out;
missing_out = out.missing_out;
w_out = out.w_out;
@rput elapsed;
@rput true_y;
@rput par_miss;
@rput missing_out;
@rput index_miss_mat;
@rput beta_out;
#@rput m_out;
#@rput c_out;
#@rput qmat_sparse_out;
@rput lambda_out;
@rput dir_plot;
@rput name;
@rput lambda;
@rput k;
@rput zmat_out;
@rput z_mat;
@rput beta_vec;
@rput y;
@rput w;
@rput beta_vec;
@rput w_out;

R"""
    library(ggplot2)
    pdf(paste(dir_plot, name, "mcmc_out_lambda.pdf", sep=""))
    hist(y)
    plot(c(y))
    abline(h = exp(beta_vec), col=2)
    plot(exp(c(w)))
    abline(h = exp(beta_vec), col=2)
    if(par_miss>0)
    {
        plot(true_y, colMeans(missing_out), pch=20)
        abline(a=0,b=1, col=2, lwd=2)
    }
    
    names_var = as.character(1:k)
    meanlambda = apply(lambda_out,c(2,3),mean)
    meanzeta = apply(zmat_out,c(2,3),mean)

    meanzeta_t = t(meanzeta)
    meanzeta_t[lower.tri(meanzeta)] = meanzeta[lower.tri(meanzeta)]
    meanzeta = t(meanzeta_t)


    hist(y)
    data_gg = data.frame(n1 = factor(rep(names_var, each = k), levels= names_var),n2 = factor(rep(names_var, times = k), levels= rev(names_var)), meanzeta = c(meanzeta), meanlambda = c(meanlambda))


    p = ggplot(data_gg, aes(x = n1, y = n2, fill = meanzeta)) +geom_tile(col ="black")+  scale_fill_gradient2(
        low = "white", 
        mid = "green", 
        high = "red", 
        midpoint = 0.5
      )+ggtitle("posterior mean of the  zeta (0 sparsa, 1 non sparsa)")

    print(p)

    
    meanzeta = z_mat

    meanzeta_t = t(meanzeta)
    meanzeta_t[lower.tri(meanzeta)] = meanzeta[lower.tri(meanzeta)]
    meanzeta = t(meanzeta_t)


    
    data_gg = data.frame(n1 = factor(rep(names_var, each = k), levels= names_var),n2 = factor(rep(names_var, times = k), levels= rev(names_var)), meanzeta = c(meanzeta), meanlambda = c(meanlambda))


    p = ggplot(data_gg, aes(x = n1, y = n2, fill = meanzeta)) +geom_tile(col ="black")+  scale_fill_gradient2(
        low = "white", 
        mid = "green", 
        high = "red", 
        midpoint = 0.5
      )+ggtitle("True  zeta (0 sparsa, 1 non sparsa)")

    print(p)

    
    p = ggplot(data_gg, aes(x = n1, y = n2, fill = meanlambda)) +geom_tile(col ="black")+ scale_fill_gradient2(
        low = "green", 
        mid = "white", 
        high = "red", 
        midpoint = .0
      )+ggtitle("posterior mean of the Inverse of the covariance matrix")

    print(p)

    meanlambda = lambda
    data_gg = data.frame(n1 = factor(rep(names_var, each = k), levels= names_var),n2 = factor(rep(names_var, times = k), levels= rev(names_var)), meanzeta = c(meanzeta), meanlambda = c(meanlambda))
    p = ggplot(data_gg, aes(x = n1, y = n2, fill = meanlambda)) +geom_tile(col ="black")+ scale_fill_gradient2(
        low = "green", 
        mid = "white", 
        high = "red", 
        midpoint = .0
      )+ggtitle("TRUE Inverse of the covariance matrix")

    print(p)

    par(mfrow=c(3,3))
    if(par_miss>0)
    {
        for(imiss in 1:dim(missing_out)[2])
        {
            plot(missing_out[,imiss], type="l", main="miss")
            abline(h=true_y[imiss], col=2, lwd=2)
        }
    }



    par(mfrow=c(3,3))
    for(j in 1:(k-1))
    {
        for(i in (j+1):k)
        {
            plot(lambda_out[,i,j], type="l", main=paste("lamda_",i,"_", j, " - TRUE=", round(lambda[i,j],3),sep=""), col= (round(lambda[i,j],7)==0)+1)
            abline(h=lambda[i,j], col=2, lwd=2)
        }
    }

    par(mfrow=c(3,3))
    for(i in 1:k)
    {
        plot(lambda_out[,i,i], type="l", main=paste("lambda_diag",i, " - TRUE=", round(lambda[i,i],3),sep=""))
        abline(h=lambda[i,i], col=2, lwd=2)
    }
    
    for(i in 1:length(beta_vec))
    {
        plot(beta_out[,i], type="l", main=paste("mu_",i, " - TRUE=", round(beta_vec[i],3),sep=""))
        abline(h=beta_vec[i], col=2, lwd=2)
    }

    dev.off()
    pdf(paste(dir_plot, name, "mcmc_out_data.pdf", sep=""))
    hist(y)
    plot(c(y))
    abline(h = exp(beta_vec), col=2)
    plot(exp(c(w)))
    abline(h = exp(beta_vec), col=2)
    par(mfrow=c(3,3))
    if(par_miss>0)
    {
        plot(true_y, colMeans(missing_out), pch=20)
        abline(a=0,b=1, col=2, lwd=2)
        
        for(iobs in 1:ncol(missing_out))
        {
            plot(missing_out[,iobs], type="l")
            abline(h = true_y[iobs], col=2, lwd= 2)
        }

    }
    par(mfrow=c(3,3))
    for(k in 1:k)
    {
        for(i in 1:dim(w_out)[3])
        {
            plot(w_out[,k,i], type="l", main = paste("w_",k,",",i))
            abline(h = w[k,i], col=2, lwd= 2)
        }
    }
    par(mfrow=c(3,3))
    for(k in 1:k)
    {
        for(i in 1:dim(w_out)[3])
        {
            plot(exp(w_out[,k,i]), type="l", main = paste("exp w_",k,",",i))
            abline(h = exp(w[k,i]), col=2, lwd= 2)
        }
    }
    dev.off()

    save.image(paste(dir_plot, name, "mcmc_out_lambda.Rdata", sep=""))
"""
