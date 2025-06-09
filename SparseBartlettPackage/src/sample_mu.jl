

function sample_mu!(iterMCMC::Int64, ydata::Matrix{Float64}, bstar_mcmc::Vector{Float64}, z_mcmc::Matrix{Float64}, model_sparse::T_SPARSE,  type_mcmc_lambda::T_HMC, psi::Matrix{Float64}, nu::Float64, r_hmc_mcmc::Vector{Float64}, q_mat_mcmc::Matrix{Float64}, q_mat_sparse_mcmc::Matrix{Float64}, lambda_mcmc::Matrix{Float64}, b_mat_mcmc::Matrix{Float64}, mu_mcmc::Matrix{Float64}, ycentered::Matrix{Float64},prior_mu::Normal{Float64}, beta_mcmc::Vector{Float64}, covariates::Array{Float64,3})  where {T_SPARSE<:TypeSparse, T_HMC<:TypeMCMC}


    k::Int64 = size(ydata,1)
    n::Int64 = size(ydata,2)
    ncov::Int64 = size(covariates,2)
    
    var_p::Matrix{Float64} = 1.0/(params(prior_mu)[2]^2) .* Matrix(I(ncov)) 
    mean_p::Vector{Float64} = params(prior_mu)[1]/(params(prior_mu)[2]^2) .* ones(Float64,ncov)
    for i = 1:size(ydata,2)
        var_p .+= transpose(covariates[:,:,i])*lambda_mcmc*covariates[:,:,i]
        mean_p .+= transpose(covariates[:,:,i])*lambda_mcmc*ydata[:,i]
    end
    #var_p = inv(cholesky(Symmetric(var_p)))
    var_p = inv(Symmetric(var_p))
    mean_p = var_p*mean_p




    beta_mcmc .= rand(MvNormal(mean_p, var_p))
    for i = 1:n
        mu_mcmc[:,i] .= covariates[:,:,i]*beta_mcmc
        ycentered[:,i] .= ydata[:,i] - mu_mcmc[:,i]
    end
    return nothing

end

#function sample_q!(iterMCMC::Int64, ydata::Matrix{Float64}, bstar_mcmc::Vector{Float64}, z_mcmc::Matrix{Int64}, model_sparse::T_SPARSE, model_psi::T_PSI, type_mcmc::HMCRJ, psi::Matrix{Float64}, nu::Float64, r_hmc_mcmc::Vector{Float64},
#    omega_mcmc::Vector{Float64},
#    gamma_mcmc::Vector{Float64},
#    q_mcmc::Vector{Float64}, q_mat_mcmc::Matrix{Float64}
#    )  where {T_SPARSE<:TypeSparse, T_PSI<:TypeS}

#    deltamax::Float64 = 1000.0
#    d_b::Vector{Float64} = zeros(Float64, size(bstar_mcmc,1))
#    log_dens = x -> compute_log_dens_hmc(x,z_mcmc, ydata, psi, nu)


#    find_reasonable_epsilon(bstar_mcmc, log_dens, d_b, type_mcmc.epsilon)

#    bstar_mcmc_minus::Vector{Float64}  = deepcopy(bstar_mcmc)
#    bstar_mcmc_plus::Vector{Float64}  = deepcopy(bstar_mcmc)
#    bstar_mcmc_start::Vector{Float64}  = deepcopy(bstar_mcmc)

#    r_hmc_mcmc_minus::Vector{Float64}  = deepcopy(r_hmc_mcmc)
#    r_hmc_mcmc_plus::Vector{Float64}  = deepcopy(r_hmc_mcmc)
#    r_hmc_mcmc_start::Vector{Float64}  = deepcopy(r_hmc_mcmc)

#    utree::Float64 = 0.0
#    jtree::Int64 = 0
#    ntree::Int64 = 1
#    stree::Int64 = 1
#    ntree_prime::Int64 = 1
#    stree_prime::Int64 = 1
#    vtree::Int64 = 1
#    for i_hmc in 1:type_mcmc.iter[1]
        
#        r_hmc_mcmc[:] = rand(Normal(0.0,1.0), size(bstar_mcmc,1))

#        utree = rand(Uniform(0.0, exp(log_dens(bstar_mcmc) -0.5 * sum(r_hmc_mcmc.^2.0) )))
#        jtree = 0
#        ntree = 1
#        stree = 1

#        bstar_mcmc_minus[:] = bstar_mcmc[:]
#        bstar_mcmc_plus[:] = bstar_mcmc[:]
#        bstar_mcmc_start[:] = bstar_mcmc[:]

#        r_hmc_mcmc_minus[:] = r_hmc_mcmc[:]
#        r_hmc_mcmc_plus[:] = r_hmc_mcmc[:]
#        r_hmc_mcmc_start[:] = r_hmc_mcmc[:]


#        while stree == 1
            
#            vtree = rand([-1, 1])

#            if vtree == -1
#                println("v-1")
#                ntree_prime, stree_prime = build_tree_hmc(log_dens, d_b , utree, vtree, jtree, type_mcmc.epsilon[1], deltamax, bstar_mcmc_minus, bstar_mcmc_minus, bstar_mcmc_minus, bstar_mcmc_start, r_hmc_mcmc_minus, r_hmc_mcmc_minus, r_hmc_mcmc_minus, r_hmc_mcmc_start)
        
#                if stree_prime == 1
#                    println("s1 " , ntree_prime,"",ntree)
#                    if rand(Uniform(0.0,1.0)) < ntree_prime/ntree
#                        println("Acc")
#                        bstar_mcmc[:] = bstar_mcmc_minus[:]
#                    end
#                end
                

#            else
#                println("v1")
#                ntree_prime, stree_prime = build_tree_hmc(log_dens, d_b , utree, vtree, jtree, type_mcmc.epsilon[1], deltamax, bstar_mcmc_plus, bstar_mcmc_plus, bstar_mcmc_plus, bstar_mcmc_start, r_hmc_mcmc_plus, r_hmc_mcmc_plus, r_hmc_mcmc_plus, r_hmc_mcmc_start)

#                if stree_prime == 1
#                    println("s1 " , ntree_prime,"",ntree)
#                    if rand(Uniform(0.0,1.0)) < ntree_prime/ntree
#                        println("Acc")
#                        bstar_mcmc[:] = bstar_mcmc_plus[:]
#                    end
#                end
                

#            end
#            ntree = ntree + ntree_prime
#            stree= stree_prime * Int64( sum((bstar_mcmc_plus .- bstar_mcmc_minus).*r_hmc_mcmc_minus) >=0.0 ) * Int64( sum((bstar_mcmc_plus .- bstar_mcmc_minus).*r_hmc_mcmc_plus)>=0.0 )
#            jtree += 1
#        end
        
#    end


#    hgamma::Int64 = 1
#    homega::Int64 = 1
#    h1::Int64 = 1
#    k::Int64 = size(q_mat_mcmc,1) 
#    for j = 1:(k-1)
        
#        q_mcmc[h1] = exp(bstar_mcmc[h1])
#        gamma_mcmc[hgamma] = q_mcmc[h1]
#        q_mat_mcmc[j,j] = q_mcmc[h1]

#        h1 += 1
#        hgamma += 1
#        for i = (j+1):k
#            q_mcmc[h1] = bstar_mcmc[h1]
#            omega_mcmc[homega] = bstar_mcmc[h1]
#            q_mat_mcmc[i,j] = q_mcmc[h1]
#            h1 += 1
#            homega += 1
#        end
#    end
#    q_mcmc[h1] = exp(bstar_mcmc[h1])
#    gamma_mcmc[hgamma] = q_mcmc[h1]
#    q_mat_mcmc[k,k] = q_mcmc[h1]


#end