
#function compute_log_dens_hmc(bstar_mcmc::Vector{Float64},z_mcmc::Matrix{Float64}, ydata::Matrix{Float64}, psi::Matrix{Float64}, nu::Float64, molt_const::Float64)::Float64

#    ret::Float64  = 0.0
#    k::Int64 = size(z_mcmc,1)
#    n::Int64 = size(ydata,2)

#    bmat::Matrix{Float64} = zeros(Float64, k,k) 
#    qmat::Matrix{Float64} = zeros(Float64, k,k) 
#    qstar::Matrix{Float64} = zeros(Float64, k,k) 


#    ###### B
#    h::Int64 = 1    
#    ## j = 1
#    bmat[1,1] = exp(bstar_mcmc[h])
#    h += 1
#    for i = 2:k
#        bmat[i,1] = bstar_mcmc[h] 
#        h += 1
#    end
#    for j = 2:(k-1)
#        bmat[j,j] = exp(bstar_mcmc[h])
#        h += 1
#        for i = (j+1):k
#            bmat[i,j] = bstar_mcmc[h]
#            h += 1
#        end
#    end
#    bmat[k,k] = exp(bstar_mcmc[h])
    
#    ###### Q = PsiB
#    for j = 1:k
#        for i = j:k
#            qmat[i,j] = sum(psi[i,j:k] .* bmat[j:k,j]) # mettere solo in funzione di bstar_mcmc
#        end
#    end    
    
#    ##### Q with zero in lambda
#    qstar[1,1] = z_mcmc[1,1]*qmat[1,1]
#    for i = 2:k
#        qstar[i,1] = z_mcmc[i,1] * qmat[i,1] 
#    end
#    for j = 2:(k-1)
#        qstar[j,j] = z_mcmc[j,j]*qmat[j,j] 
#        for i = (j+1):k
#            qstar[i,j] = z_mcmc[i,j] * qmat[i,j] + (1.0-z_mcmc[i,j]) * (- sum(qstar[i,1:(j-1)] .* qstar[j,1:(j-1)]) / qstar[j,j])
#        end
#    end
#    qstar[k,k] = z_mcmc[k,k]*qmat[k,k]

#    ## like
#    ret += -molt_const
#    for iobs = 1:n
#        for j = 1:k
#            ret += log(qstar[j,j]) -0.5*( sum(ydata[j:k,iobs] .* qstar[j:k,j] ))^2.0
#        end
#    end

#    ## prior
#    for j = 1:(k-1)
#        for i = (j+1):k
#            ret += -0.5 * bmat[i,j]^2
#        end
#    end
#    for j = 1:k
#        ret +=  2.0*( (nu-j+1)/2.0 - 1.0)*log(bmat[j,j]) -bmat[j,j]^2/2.0 + 2.0*log(bmat[j,j]) 
#    end
    
#    return ret

#end
#b = [-26989.112490382933, -25986.774280068475, -23948.06290872419, -197738.69332942477, -39078.01956592724, -30053.323895443213, -1.0616057474626538e6, -86495.72778156682, -918050.7584891304, -176894.90826880085, -143405.67614218537, -393783.8572413009, -1.066550494777848e6, -202745.7684574463, -165385.70300710527, -2.511054696096925e7, -1.6650182952587055e6, -1.464197070335808e6, -872892.2546959647, -256567.7486609355, -737570.0813001486]

#a = [-26989.112490382933, -25986.774280068475, -23948.06290872419, -197738.69332942477, -39078.01956592724, -30053.323895443213, -1.0616057474626538e6, -86495.72778156682, -918050.7584891304, -176894.90826880085, -143405.67614218537, -393783.8572413009, -1.066550494777848e6, -202745.7684574463, -165385.70300710527, -2.511054696096925e7, -1.6650182952587055e6, -1.464197070335808e6, -872892.2546959647, -256567.7486609355, -737570.0813001486]

#-5.3890560989306495
#ERROR:-328375.1755994872
function compute_log_dens_hmc(bstar_mcmc::Vector{Float64}, z_mcmc::Matrix{Float64}, ycentered::Matrix{Float64}, psi::Matrix{Float64}, nu::Float64, molt_const::Float64, n_uno_vec::Vector{Float64} )::Float64

    ret::Float64  = 0.0
    k::Int64 = size(z_mcmc,1)
    n::Int64 = size(ycentered,2)

    #bmat::Matrix{Float64} = zeros(Float64, k,k) 
    #qmat::Matrix{Float64} = zeros(Float64, k,k) 
     qstar::Matrix{Float64} = zeros(Float64, k,k) 
 
    
    
    ret += -molt_const
    ###### B
    h::Int64 = 1  
    
    ## j = 1
    qstar[1,1] = exp(bstar_mcmc[h])
    
    h += 1
    
    for i = 2:k
        #n_uno += (1-z_mcmc[i,1])
        qstar[i,1] = bstar_mcmc[h] 
        ret += -0.5 * bstar_mcmc[h]^2.0
        h += 1
    end
    ret +=  2.0*( (nu-n_uno_vec[1])/2.0 - 1.0)*log(qstar[1,1]) -qstar[1,1]^2/2.0 + 2.0*log(qstar[1,1]) 
    for j = 2:(k-1)
        qstar[j,j] = exp(bstar_mcmc[h])
        
        #for i = 1:(j-1)
        #    n_uno += (1-z_mcmc[j,i])
        #end
        
        h += 1
        for i = (j+1):k
            #n_uno += (1-z_mcmc[i,j])
            qstar[i,j] = bstar_mcmc[h]
            ret += -0.5 * bstar_mcmc[h]^2.0
            h += 1
        end
        ret += 2.0 * ((nu - n_uno_vec[j]) / 2.0 - 1.0) * log(qstar[j, j]) - qstar[j, j]^2 / 2.0 + 2.0 * log(qstar[j, j])
    end
    qstar[k,k] = exp(bstar_mcmc[h])
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

    ## like
    
    for iobs = 1:n
        for j = 1:k
            ret += log(qstar[j,j]) -0.5*( sum(ycentered[j:k,iobs] .* qstar[j:k,j] ))^2.0
        end
    end

    
    
    return ret

end


function sample_q!(iterMCMC::Int64, ydata::Matrix{Float64}, bstar_mcmc::Vector{Float64}, z_mcmc::Matrix{Float64}, model_sparse::T_SPARSE,  type_mcmc_lambda::T_HMC, psi::Matrix{Float64}, nu_orig::Float64, r_hmc_mcmc::Vector{Float64}, q_mat_mcmc::Matrix{Float64}, q_mat_sparse_mcmc::Matrix{Float64}, lambda_mcmc::Matrix{Float64}, b_mat_mcmc::Matrix{Float64}, mu_mcmc::Matrix{Float64}, ycentered::Matrix{Float64})  where {T_SPARSE<:TypeSparse, T_HMC<:TypeMCMC}

    k::Int64 = size(z_mcmc, 1)
         
    n_uno = zeros(Float64,k)
    n_zero = zeros(Float64,k)
    
    # ! numeri uno per riga
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
    d_b::Vector{Float64} = zeros(Float64, size(bstar_mcmc,1))
    
    molt_const = compute_log_dens_hmc(bstar_mcmc, z_mcmc, ycentered, psi, nu, 0.0, n_uno )
    
    #if  isnan(molt_const)

    #    k::Int64 = size(z_mcmc,1)
    #    n::Int64 = size(ycentered,2)
    #    nm::Int64 = size(bstar_mcmc,1)- k
        
    #    bmat::Matrix{Float64} = zeros(Float64, k,k) 
    #    from_bstar_to_bmat(bstar_mcmc, bmat, k, n)
    #    qmat::Matrix{Float64} = zeros(Float64, k,k) 
    #    from_bmat_to_qmat(bmat, qmat,k, n, psi)
    #    qstar::Matrix{Float64} = zeros(Float64, k,k)
    #    from_qmat_to_qstar_init(qmat, qstar, k, n, psi, z_mcmc)


    #    like_acc::Vector{Float64} = zeros(Float64,n)
    #    for j = 1:k
    #        compute_like_for_z(j, n, qstar, like_acc, z_mcmc, k ,ycentered)
    #    end
    #    println(like_acc)
    #    println("a  = ", qmat)
    #end
    log_dens = x -> compute_log_dens_hmc(x, z_mcmc, ycentered, psi, nu, molt_const, n_uno )

    #println(molt_const)
    #error("")
    gradient!(Reverse, d_b, log_dens, bstar_mcmc)
    #println(log_dens)
    #println(sum(d_b))
    #error("-328375.1755994872")
    
    add_const::Float64 = 0.5*size(bstar_mcmc,1)
    #add_const = 0.0
    hmc_sample(log_dens, d_b, bstar_mcmc, r_hmc_mcmc, type_mcmc_lambda, iterMCMC,add_const )
    
    #dual_av_no_uturn(log_dens, d_b,   type_mcmc_lambda.epsilon, 1000.0, bstar_mcmc,  r_hmc_mcmc, 1, 0.9, 10)
    
    from_bstar_to_b_mat(bstar_mcmc, b_mat_mcmc)
    from_bmat_to_all_q(b_mat_mcmc,  q_mat_mcmc, q_mat_sparse_mcmc, z_mcmc,psi)
    lambda_mcmc[:,:] = q_mat_sparse_mcmc*transpose(q_mat_sparse_mcmc)
    
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