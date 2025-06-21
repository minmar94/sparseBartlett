

dir = "/Users/gianlucamastrantonio/Library/CloudStorage/GoogleDrive-mastrantonio.gluca@gmail.com/Shared drives/Size and Shape/dataset/"
dir_plot = "/Users/gianlucamastrantonio/Library/CloudStorage/GoogleDrive-mastrantonio.gluca@gmail.com/Shared drives/Size and Shape/dataset/results/"
# 3 e 15
name = "doubs_indipendent_species"
using Pkg
Pkg.activate(dir)


using MKL
using RCall
using Distributions
using PDMats
using LinearAlgebra
using ToggleableAsserts
using Enzyme
using Random
using SparseIW
Random.seed!(1);
## sim
R"""

dataset = read.csv("/Users/gianlucamastrantonio/Library/CloudStorage/GoogleDrive-mastrantonio.gluca@gmail.com/Shared drives/Size and Shape/dataset/doubs.csv", header = F)
summary(dataset)
#dataset = rowSums(dataset)

site = dataset[1,]
species = dataset[,1]

dataset = dataset[-1,]
dataset = dataset[,-1]
dataset = as.matrix(dataset)

dataset <- t(dataset[,-8])
"""

@rget dataset;
dataset = 1.0 .* dataset
n::Int64 = size(dataset, 2)
k::Int64 = size(dataset, 1)


##### ##### ##### ##### ##### ##### 
##### MCMC 
##### ##### ##### ##### ##### ##### 
z_mat_init = ones(Float64, k, k)

molt = 1

nm = Int64((k * k - k) / 2)
#mcmc_par = (iter = 30, burnin = 10, thin = 1)
mcmc_par = (iter=10000 * molt, burnin=8000 * molt, thin=1 * molt)

#mcmc_par = (iter=10000 * molt, burnin=8000 * molt, thin=1 * molt)
covariate = ones(Float64, k, 1, n)

#mcmc_par = (iter = 100 , burnin = 80 , thin = 1  )
#elapsed = @elapsed begin
    out = mcmc_glm(
    ydata=dataset,
        iterations=mcmc_par,
    covariates=covariate,
        model_sparse=GenSparse(Int64.(nm:(-1):1) ./ sum(Int64.(nm:(-1):1))),
        model_psi=GeneralS(),
        #type_mcmc_lambda = BaseHMC(;epsilon = 1.0, L = 10, M = 10, restart = true),
        #type_mcmc_lambda = HMC_DualAv(;epsilon = 1.0, L = 10, M = 10, restart = false, delta = 0.5, lambda = 0.001, iter_adapt = 10000000   ),
        type_mcmc_lambda=HMC_DualAv_NoUTurn(; epsilon=10.000, L=10, M=10, restart=false, delta=0.5, iter_adapt=10000000),
        #type_mcmc = HMC(0.8),

        prior_lambda=Wishart(1.0 * (k + 3 - 1), 1.0 .* Matrix(I(k))),
        #prior_lambda = Wishart(1.0*(k+3-1), (1.0/(k+3-1)).*Matrix(I(k))),
        #prior_lambda = Wishart(k-0.5, Matrix(I(k))./(k-0.5)),
        prior_mu=Normal(0.0, 1000.0^0.5), c_init=ones(Float64, k),
        m_init=zeros(Float64, nm),
    beta_init=zeros(Float64, size(covariate, 2)),
        prob_zeta=ones(Float64, k, k) * 0.5,
        #m_init = vec_m,
        #c_init = vec_c,#[sqrt(rand(Chi((k+2)-i+1))) for i = 1:k ],
        z_init=Float64.(z_mat_init),
        w_init=zeros(Float64, k, n),
        like=Poisson(1.0),
        iter_start_hmc=2,
        iter_start_zeta=5,
        
        #missing_index=index_miss_mat
    )
#end
m_out = out.m_out;
c_out = out.c_out;
qmat_sparse_out = out.q_mat_sparse_out;
lambda_out = out.lambda_out;
zmat_out = out.zmat_out;
beta_out = out.beta_out;
y = dataset;


@rput y;
#@rput par_miss;
#@rput missing_out;
#@rput index_miss_mat;
@rput beta_out;
#@rput m_out;
#@rput c_out;
#@rput qmat_sparse_out;
@rput lambda_out;
@rput dir_plot;
@rput name;
#@rput lambda;
@rput k;
@rput zmat_out;




R"""
    library(ggplot2)
    pdf(paste(dir_plot, name, "mcmc_out_lambda.pdf", sep=""))
    #hist(y)
    #plot(c(y))
    #if(par_miss>0)
    #{
    #    plot(true_y, colMeans(missing_out), pch=20)
    #    abline(a=0,b=1, col=2, lwd=2)
    #}
    
    names_var = as.character(1:k)
    meanlambda = apply(lambda_out,c(2,3),mean)
    meanzeta = apply(zmat_out,c(2,3),mean)

    meanzeta_t = t(meanzeta)
    meanzeta_t[lower.tri(meanzeta)] = meanzeta[lower.tri(meanzeta)]
    meanzeta = t(meanzeta_t)


    #hist(y)
    data_gg = data.frame(n1 = factor(rep(names_var, each = k), levels= names_var),n2 = factor(rep(names_var, times = k), levels= rev(names_var)), meanzeta = c(meanzeta), meanlambda = c(meanlambda))


    p = ggplot(data_gg, aes(x = n1, y = n2, fill = meanzeta)) +geom_tile(col ="black")+  scale_fill_gradient2(
        low = "white", 
        mid = "green", 
        high = "red", 
        midpoint = 0.5
      )+ggtitle("posterior mean of the  zeta (0 sparsa, 1 non sparsa)")

    print(p)

    
    #meanzeta = z_mat

    #meanzeta_t = t(meanzeta)
    #meanzeta_t[lower.tri(meanzeta)] = meanzeta[lower.tri(meanzeta)]
    #meanzeta = t(meanzeta_t)


    
    #data_gg = data.frame(n1 = factor(rep(names_var, each = k), levels= names_var),n2 = factor(rep(names_var, times = k), levels= rev(names_var)), meanzeta = c(meanzeta), meanlambda = c(meanlambda))


    #p = ggplot(data_gg, aes(x = n1, y = n2, fill = meanzeta)) +geom_tile(col ="black")+  scale_fill_gradient2(
    #    low = "white", 
    #    mid = "green", 
    #    high = "red", 
    #    midpoint = 0.5
    #  )+ggtitle("True  zeta (0 sparsa, 1 non sparsa)")

    #print(p)

    
    p = ggplot(data_gg, aes(x = n1, y = n2, fill = meanlambda)) +geom_tile(col ="black")+ scale_fill_gradient2(
        low = "green", 
        mid = "white", 
        high = "red", 
        midpoint = .0
      )+ggtitle("posterior mean of the Inverse of the covariance matrix")

    print(p)

    #meanlambda = lambda
    #data_gg = data.frame(n1 = factor(rep(names_var, each = k), levels= names_var),n2 = factor(rep(names_var, times = k), levels= rev(names_var)), meanzeta = c(meanzeta), meanlambda = c(meanlambda))
    #p = ggplot(data_gg, aes(x = n1, y = n2, fill = meanlambda)) +geom_tile(col ="black")+ scale_fill_gradient2(
    #    low = "green", 
    #    mid = "white", 
    #    high = "red", 
    #    midpoint = .0
    #  )+ggtitle("TRUE Inverse of the covariance matrix")

    #print(p)

    #par(mfrow=c(3,3))
    #if(par_miss>0)
    #{
    #    for(imiss in 1:dim(missing_out)[2])
    #    {
    #        plot(missing_out[,imiss], type="l", main="miss")
    #        abline(h=true_y[imiss], col=2, lwd=2)
    #    }
    #}



    par(mfrow=c(3,3))
    for(j in 1:(k-1))
    {
        for(i in (j+1):k)
        {
            plot(lambda_out[,i,j], type="l", main=paste("lamda_",i,"_", j, sep=""))
            #abline(h=lambda[i,j], col=2, lwd=2)
        }
    }

    par(mfrow=c(3,3))
    for(i in 1:k)
    {
        plot(lambda_out[,i,i], type="l", main=paste("lambda_diag",i,sep=""))
        #abline(h=lambda[i,i], col=2, lwd=2)
    }
    
    for(i in 1:dim(beta_out)[2])
    {
        plot(beta_out[,i], type="l", main=paste("mu_",i, sep=""))
        #abline(h=beta_vec[i], col=2, lwd=2)
    }

    dev.off()
    #pdf(paste(dir_plot, name, "mcmc_out_data.pdf", sep=""))
    #hist(y)
    #plot(c(y))
    #par(mfrow=c(3,3))
    #if(par_miss>0)
    #{
    #    plot(true_y, colMeans(missing_out), pch=20)
    #    abline(a=0,b=1, col=2, lwd=2)
        
    #    for(iobs in 1:ncol(missing_out))
    #    {
    #        plot(missing_out[,iobs], type="l")
    #        abline(h = true_y[iobs], col=2, lwd= 2)
    #    }

    #}
    #dev.off()

    save.image(paste(dir_plot, name, "mcmc_out_lambda.Rdata", sep=""))
"""
