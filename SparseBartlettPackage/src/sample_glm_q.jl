
function compute_log_dens_hmc(bwstar_mcmc::Vector{Float64}, z_mcmc::Matrix{Float64}, ydata::Matrix{Float64}, psi::Matrix{Float64}, nu::Float64, w_mcmc::Matrix{Float64}, mu_mcmc::Matrix{Float64}, like::Bernoulli, molt_const::Float64, n_uno_vec::Vector{Float64})::Float64

    ret::Float64  = 0.0
    k::Int64 = size(z_mcmc,1)
    n::Int64 = size(ydata,2)
    nm::Int64 = (k^2-k)/2+k
    #bmat::Matrix{Float64} = zeros(Float64, k,k) 
    #qmat::Matrix{Float64} = zeros(Float64, k,k) 
    qstar::Matrix{Float64} = zeros(Float64, k,k) 
    wmcmc::Matrix{Float64} = zeros(Float64, k,n) 

    

    ret += -molt_const
    ###### B
    h::Int64 = 1    
    ## j = 1
    qstar[1,1] = exp(bwstar_mcmc[h])
    ret += 2.0 * ((nu - n_uno_vec[1]) / 2.0 - 1.0) * log(qstar[1, 1]) - qstar[1, 1]^2 / 2.0 + 2.0 * log(qstar[1, 1])
    h += 1
    for i = 2:k
        qstar[i,1] = bwstar_mcmc[h] 
        ret += -0.5 * bwstar_mcmc[h]^2.0
        h += 1
    end
    for j = 2:(k-1)
        qstar[j,j] = exp(bwstar_mcmc[h])
        ret += 2.0 * ((nu - n_uno_vec[j]) / 2.0 - 1.0) * log(qstar[j, j]) - qstar[j, j]^2 / 2.0 + 2.0 * log(qstar[j, j])
        h += 1
        for i = (j+1):k
            qstar[i,j] = bwstar_mcmc[h]
            ret += -0.5 * bwstar_mcmc[h]^2.0
            h += 1
        end
    end
    qstar[k,k] = exp(bwstar_mcmc[h])
    ret += 2.0 * ((nu - n_uno_vec[k]) / 2.0 - 1.0) * log(qstar[k, k]) - qstar[k, k]^2 / 2.0 + 2.0 * log(qstar[k, k])
    
    
    ###### Q = PsiB
    for j = 1:k
        for i = j:k
            qstar[i,j] = sum(psi[i,j:k] .* qstar[j:k,j]) # mettere solo in funzione di bstar_mcmc
        end
    end    
    
    ##### Q with zero in lambda
    qstar[1,1] = z_mcmc[1,1]*qstar[1,1]
    for i = 2:k
        qstar[i,1] = z_mcmc[i,1] * qstar[i,1] 
    end
    for j = 2:(k-1)
        qstar[j,j] = z_mcmc[j,j]*qstar[j,j] 
        for i = (j+1):k
            qstar[i,j] = z_mcmc[i,j] * qstar[i,j] + (1.0-z_mcmc[i,j]) * (- sum(qstar[i,1:(j-1)] .* qstar[j,1:(j-1)]) / qstar[j,j])
        end
    end
    qstar[k,k] = z_mcmc[k,k]*qstar[k,k]

    for iobs = 1:n
        for j = 1:k
            wmcmc[j,iobs] = bwstar_mcmc[nm +  (iobs-1)*k + j]
        end
    end

    
    ## like
    
    for iobs = 1:n
        for j = 1:k
            ret += log(qstar[j,j]) -0.5*( sum((wmcmc[j:k,iobs]-mu_mcmc[j:k, iobs]) .* qstar[j:k,j] ))^2.0
        end
    end
    for iobs = 1:n
        for j = 1:k
            ret += ydata[j,iobs]*wmcmc[j,iobs] - log(1.0+exp(wmcmc[j,iobs]))
        end
    end

    
    
    return ret

end


function compute_log_dens_hmc(bwstar_mcmc::Vector{Float64}, z_mcmc::Matrix{Float64}, ydata::Matrix{Float64}, psi::Matrix{Float64}, nu::Float64, w_mcmc::Matrix{Float64}, mu_mcmc::Matrix{Float64}, like::Poisson, molt_const::Float64, n_uno_vec::Vector{Float64})::Float64

    ret::Float64  = 0.0
    k::Int64 = size(z_mcmc,1)
    n::Int64 = size(ydata,2)
    nm::Int64 = (k^2-k)/2+k
    #bmat::Matrix{Float64} = zeros(Float64, k,k) 
    #qmat::Matrix{Float64} = zeros(Float64, k,k) 
    qstar::Matrix{Float64} = zeros(Float64, k,k) 
    wmcmc::Matrix{Float64} = zeros(Float64, k,n) 

    

    ret += -molt_const
    ###### B
    h::Int64 = 1    
    ## j = 1
    qstar[1,1] = exp(bwstar_mcmc[h])
    ret += 2.0 * ((nu - n_uno_vec[1]) / 2.0 - 1.0) * log(qstar[1, 1]) - qstar[1, 1]^2 / 2.0 + 2.0 * log(qstar[1, 1])
    h += 1
    for i = 2:k
        qstar[i,1] = bwstar_mcmc[h] 
        ret += -0.5 * bwstar_mcmc[h]^2.0
        h += 1
    end
    for j = 2:(k-1)
        qstar[j,j] = exp(bwstar_mcmc[h])
        ret += 2.0 * ((nu - n_uno_vec[j] ) / 2.0 - 1.0) * log(qstar[j, j]) - qstar[j, j]^2 / 2.0 + 2.0 * log(qstar[j, j])
        h += 1
        for i = (j+1):k
            qstar[i,j] = bwstar_mcmc[h]
            ret += -0.5 * bwstar_mcmc[h]^2.0
            h += 1
        end
    end
    qstar[k,k] = exp(bwstar_mcmc[h])
    ret +=  2.0*( (nu - n_uno_vec[k]) / 2.0 - 1.0)*log(qstar[k,k]) -qstar[k,k]^2/2.0 + 2.0*log(qstar[k,k]) 
    
    
    ###### Q = PsiB
    for j = 1:k
        for i = j:k
            qstar[i,j] = sum(psi[i,j:k] .* qstar[j:k,j]) # mettere solo in funzione di bstar_mcmc
        end
    end    
    
    ##### Q with zero in lambda
    qstar[1,1] = z_mcmc[1,1]*qstar[1,1]
    for i = 2:k
        qstar[i,1] = z_mcmc[i,1] * qstar[i,1] 
    end
    for j = 2:(k-1)
        qstar[j,j] = z_mcmc[j,j]*qstar[j,j] 
        for i = (j+1):k
            qstar[i,j] = z_mcmc[i,j] * qstar[i,j] + (1.0-z_mcmc[i,j]) * (- sum(qstar[i,1:(j-1)] .* qstar[j,1:(j-1)]) / qstar[j,j])
        end
    end
    qstar[k,k] = z_mcmc[k,k]*qstar[k,k]

    for iobs = 1:n
        for j = 1:k
            wmcmc[j,iobs] = bwstar_mcmc[nm +  (iobs-1)*k + j]
        end
    end

    
    ## like
    
    for iobs = 1:n
        for j = 1:k
            ret += log(qstar[j,j]) -0.5*( sum((wmcmc[j:k,iobs]-mu_mcmc[j:k, iobs]) .* qstar[j:k,j] ))^2.0
        end
    end
    for iobs = 1:n
        for j = 1:k
            ret += ydata[j,iobs]*wmcmc[j,iobs] - exp(wmcmc[j,iobs])
        end
    end

    #lambda^x e^-lambda
#@x log(lambda)  -log(lambda)
    
    
    return ret

end


function sample_glm_q!(iterMCMC::Int64, ydata::Matrix{Float64}, bwstar_mcmc::Vector{Float64}, z_mcmc::Matrix{Float64}, model_sparse::T_SPARSE,  type_mcmc_lambda::T_HMC, psi::Matrix{Float64}, nu_orig::Float64, r_hmc_mcmc::Vector{Float64}, q_mat_mcmc::Matrix{Float64}, q_mat_sparse_mcmc::Matrix{Float64}, lambda_mcmc::Matrix{Float64}, b_mat_mcmc::Matrix{Float64}, mu_mcmc::Matrix{Float64}, w_mcmc::Matrix{Float64}, like::DD, wcentered::Matrix{Float64})  where {T_SPARSE<:TypeSparse, T_HMC<:TypeMCMC, DD<: UnivariateDistribution}


    
    k::Int64 = size(z_mcmc, 1)

    n_uno = zeros(Float64,k)
    n_zero = zeros(Float64,k)

    # ! numeri 1 per riga
    #for j = 2:k
    #    for i = 1:(j-1)
    #        n_uno[j] += Float64(z_mcmc[j,i])
    #        n_zero[j] += Float64(1.0 - z_mcmc[j,i])
    #    end
    #end
    #nu::Float64 = nu_orig -  k+1 + maximum(n_uno)

    # ! numeri 1 per colonna
    for i = 1:(k-1)
        for j = (i+1):(k)
            n_uno[i] -= Float64(z_mcmc[j,i])
            #n_zero[j] += Float64(1.0 - z_mcmc[j,i])
        end
    end
    nu::Float64 = nu_orig -  k+1
    #deltamax::Float64 = 1000.0
    d_b::Vector{Float64} = zeros(Float64, size(bwstar_mcmc,1))

    molt_const = compute_log_dens_hmc(bwstar_mcmc, z_mcmc, ydata, psi, nu, w_mcmc, mu_mcmc, like, 0.0, n_uno)
    
    log_dens = x -> compute_log_dens_hmc(x, z_mcmc, ydata, psi, nu, w_mcmc, mu_mcmc, like, molt_const, n_uno)

    gradient!(Reverse, d_b, log_dens, bwstar_mcmc)
    
    #println(d_b)
    #error("R")
    
    add_const::Float64 = 0.5*size(bwstar_mcmc,1)
    #add_const = 0.0
    hmc_sample(log_dens, d_b, bwstar_mcmc, r_hmc_mcmc, type_mcmc_lambda, iterMCMC,add_const )
    
    #dual_av_no_uturn(log_dens, d_b,   type_mcmc_lambda.epsilon, 1000.0, bstar_mcmc,  r_hmc_mcmc, 1, 0.9, 10)
    
    from_bstar_to_b_mat(bwstar_mcmc, b_mat_mcmc)
    from_bmat_to_all_q(b_mat_mcmc,  q_mat_mcmc, q_mat_sparse_mcmc, z_mcmc,psi)
    lambda_mcmc[:,:] = q_mat_sparse_mcmc*transpose(q_mat_sparse_mcmc)
    #k::Int64 = size(z_mcmc,1)
    n::Int64 = size(ydata,2)
    nm::Int64 = (k^2-k)/2+k
    for iobs = 1:n
        for j = 1:k
            w_mcmc[j,iobs] = bwstar_mcmc[nm +  (iobs-1)*k + j ]
            wcentered[j,iobs] = w_mcmc[j,iobs] -mu_mcmc[j,iobs]
        end
    end


    return nothing

end
