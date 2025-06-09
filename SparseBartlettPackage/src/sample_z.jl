function from_bstar_to_bmat(bstar_mcmc::Vector{Float64}, bmat::Matrix{Float64}, k::Int64, n::Int64)

    h::Int64 = 1    
    ## j = 1
    bmat[1,1] = exp(bstar_mcmc[h])
    h += 1
    for i = 2:k
        bmat[i,1] = bstar_mcmc[h] 
        h += 1
    end
    for j = 2:(k-1)
        bmat[j,j] = exp(bstar_mcmc[h])
        h += 1
        for i = (j+1):k
            bmat[i,j] = bstar_mcmc[h]
            h += 1
        end
    end
    bmat[k,k] = exp(bstar_mcmc[h])

    return nothing


end

function from_bmat_to_qmat(bmat::Matrix{Float64}, qmat::Matrix{Float64},k::Int64, n::Int64, psi::Matrix{Float64})

    for j = 1:k
        for i = j:k
            qmat[i,j] = sum(psi[i,j:k] .* bmat[j:k,j]) # mettere solo in funzione di bstar_mcmc
        end
    end  
    
    return nothing


end


function from_qmat_to_qstar_init(qmat::Matrix{Float64}, qstar::Matrix{Float64},k::Int64, n::Int64, psi::Matrix{Float64}, z_mcmc::Matrix{Float64})


    
    qstar[1,1] = z_mcmc[1,1]*qmat[1,1]
    for i = 2:k
        from_qmat_to_qstar(i, 1, qmat, qstar,k, n, psi, z_mcmc)
    end
    for j = 2:(k-1)
        qstar[j,j] = z_mcmc[j,j]*qmat[j,j] 
        for i = (j+1):k
            from_qmat_to_qstar(i, j, qmat, qstar,k, n, psi, z_mcmc)
        end
    end
    qstar[k,k] = z_mcmc[k,k]*qmat[k,k]

    
    return nothing


end


function from_qmat_to_qstar_triangularindex(i::Int64, j::Int64, qmat::Matrix{Float64}, qstar::Matrix{Float64},k::Int64, n::Int64, psi::Matrix{Float64}, z_mcmc::Matrix{Float64})


    
    #for irow = (j+1):k
    #    for jcol = j:(k-1)
    #        from_qmat_to_qstar(irow, jcol, qmat, qstar,k, n, psi, z_mcmc)
    #    end
    #end

    from_qmat_to_qstar(i, j, qmat, qstar,k, n, psi, z_mcmc)
    if (j+1) < i
        for jcol  = (j+1):(i-1)
            for irow = i:k
                from_qmat_to_qstar(irow, jcol, qmat, qstar,k, n, psi, z_mcmc)
            end
        end
    end
    if i<k
        for jcol = i:(k-1)                            
            for irow = (jcol+1):k
                from_qmat_to_qstar(irow, jcol, qmat, qstar,k, n, psi, z_mcmc)
            end
        end
    end
    
    return nothing


end

function from_qmat_to_qstar(i::Int64, j::Int64, qmat::Matrix{Float64}, qstar::Matrix{Float64},k::Int64, n::Int64, psi::Matrix{Float64}, z_mcmc::Matrix{Float64})


    
    if j == 1
        qstar[i,1] = z_mcmc[i,1] * qmat[i,1] 
    else
        if z_mcmc[i,j] == 1.0
            qstar[i,j] = qmat[i,j]
        else
            qstar[i,j] = 0.0
            for jcol = 1:(j-1)
                qstar[i,j] += qstar[i,jcol] .* qstar[j,jcol]
            end
            qstar[i,j] = - qstar[i,j]/qstar[j,j]
            #qstar[i,j] = - sum(qstar[i,1:(j-1)] .* qstar[j,1:(j-1)]) / qstar[j,j]
        end
        #qstar[i,j] = z_mcmc[i,j] * qmat[i,j] + (1.0-z_mcmc[i,j]) * (- sum(qstar[i,1:(j-1)] .* qstar[j,1:(j-1)]) / qstar[j,j])
    end
    #qstar[k,k] = z_mcmc[k,k]*qmat[k,k]

    
    return nothing


end

function compute_like_for_z(j::Int64, n::Int64,qstar::Matrix{Float64}, like::Vector{Float64}, z_mcmc::Matrix{Float64} ,k::Int64, ydata::Matrix{Float64})

    like[j] = 0.0
    for iobs = 1:n
        like[j] += compute_like_for_z_colj(iobs, j, n, qstar, like , z_mcmc, k, ydata )
    end

end

function compute_like_for_z_colj(iobs::Int64,  j::Int64, n::Int64, qstar::Matrix{Float64}, like_acc::Vector{Float64} , z_mcmc::Matrix{Float64} ,k::Int64, ydata::Matrix{Float64})::Float64

    ret::Float64 = 0.0
    for jcol = j:k
        ret += ydata[jcol,iobs] .* qstar[jcol,j]
    end
    return log(qstar[j,j]) -0.5*ret^2.0
    #return log(qstar[j,j]) -0.5*( sum(ydata[j:k,iobs] .* qstar[j:k,j] ))^2.0
end

function sample_zeta!(iterMCMC::Int64, ydata::Matrix{Float64}, bstar_mcmc::Vector{Float64}, z_mcmc::Matrix{Float64}, model_sparse::T_SPARSE,  type_mcmc_lambda::T_HMC, psi::Matrix{Float64}, nu::Float64, r_hmc_mcmc::Vector{Float64}, q_mat_mcmc::Matrix{Float64}, q_mat_sparse_mcmc::Matrix{Float64}, lambda_mcmc::Matrix{Float64}, b_mat_mcmc::Matrix{Float64},ycentered::Matrix{Float64},prob_zeta_mcmc::Matrix{Float64})  where {T_SPARSE<:TypeSparse, T_HMC<:TypeMCMC}

    sum_zero::Int64 = abs(sum(z_mcmc  .- 1)) +  1
    
    
    #n_sample::Int64 = size(model_sparse.prob_vec_zeros)[1]
    k::Int64 = size(z_mcmc,1)
    n::Int64 = size(ycentered,2)
    nm::Int64 = size(bstar_mcmc,1)- k
    hc::Int64 = 0
    #hc_vec = zeros(Int64,k)
    #hh::Int64 = 1
    #for j = 1:k
    #    hc_vec[j] = hh
    #    hh = hh + (k-j) + 1
    #end
    #println("start")
    #println([nu, Float64(k)])
    #n_zero_vec = zeros(Float64, k)
    #for j = 1:(k-1)
    #    for i = (j+1):k
    #        n_zero_vec[j] += Float64((1.0 - z_mcmc[i, j]))
    #    end
    #    #@assert     n_zero_vec[j] >= 0
    #    #@assert     n_zero_vec[j] <= k
    #end
    n_uno_vec = zeros(Float64,k)
    n_zero_vec = zeros(Float64,k)
    # ! numeri 1 per riga
    #for j = 2:k
    #    for i = 1:(j-1)
    #        n_uno_vec[j] += Float64(z_mcmc[j,i])
    #        n_zero_vec[j] += Float64(1.0 - z_mcmc[j,i])
    #    end
    #end
    # ! numeri 1 per colonna
    for i = 1:(k-1)
        for j = (i+1):(k)
            n_uno_vec[i] += Float64(z_mcmc[j,i])
            #n_zero[j] += Float64(1.0 - z_mcmc[j,i])
        end
    end

    n_uno_vec_prop = deepcopy(n_uno_vec)
    n_zero_vec_prop = deepcopy(n_zero_vec)
    #println(n_zero_vec)

    samp::Int64 = 0
    z_prop::Matrix{Float64} = deepcopy(z_mcmc)

    bmat::Matrix{Float64} = zeros(Float64, k,k) 
    from_bstar_to_bmat(bstar_mcmc, bmat, k, n)
    qmat::Matrix{Float64} = zeros(Float64, k,k) 
    from_bmat_to_qmat(bmat, qmat,k, n, psi)
    qstar::Matrix{Float64} = zeros(Float64, k,k)
    from_qmat_to_qstar_init(qmat, qstar, k, n, psi, z_mcmc)
    qstar_prop = deepcopy(qstar)


    like_acc::Vector{Float64} = zeros(Float64,k)
    for j = 1:k
        compute_like_for_z(j, n, qstar, like_acc, z_mcmc, k ,ycentered)
    end
    like_prop::Vector{Float64} = deepcopy(like_acc)


    log_prob_0_1::Vector{Float64} = zeros(Float64, 2)
    nzero_0_1::Vector{Int64} = zeros(Int64, 2)

    #sum_prob_over1::Float64 = 0.0
    #for ii = 1:(k-1)
    #    for jj = (ii+1):k
    #        if z_mcmc[jj,ii] == 1
    #            #println(prob_zeta_mcmc[jj,ii])
    #            sum_prob_over1 += prob_zeta_mcmc[jj,ii]
    #        end
            
    #    end
    #end
#    println(prob_zeta_mcmc)
#    println(z_mcmc)
#    println(sum_prob_over1)
#error("")
    par0::Float64 = 0.00
    par1::Float64 = 0.00
    for s = 1:nm

        i = sample(1:k,1)[1]
        j = sample(1:k,1)[1]

        #@assert (z_mcmc[i,j] == 0) | (z_mcmc[i,j]  == 1)
        if i != j
 
            j,i = min(i,j), max(i,j)
            if z_mcmc[i,j] == 0.0
                
                # ! colonna
                #n_uno_vec_prop[i] = n_uno_vec[i] + 1.0
                #n_zero_vec_prop[i] = n_zero_vec[i] - 1.0
                # prob di mettere uno zero
                #par0 = (nu - Float64(k) + 1.0  + maximum(n_uno_vec) -n_uno_vec[i] )
                #par1 = (nu - Float64(k) + 1.0  + maximum(n_uno_vec_prop) -n_uno_vec_prop[i] )
                
                # ! riga
                n_uno_vec_prop[j] = n_uno_vec[j] + 1.0
                n_zero_vec_prop[j] = n_zero_vec[j] - 1.0
                # prob di mettere uno zero
                par0 = (nu - Float64(k) + 1.0  + n_uno_vec[j] )
                par1 = (nu - Float64(k) + 1.0  + n_uno_vec_prop[j] )
            
                
                

                log_prob_0_1[1] = log(1.0 - prob_zeta_mcmc[i, j]) + 2.0 * (par0 / 2.0 - 1.0) * log(bmat[j,j]) + (par0/ 2.0)*log(0.5) - loggamma(par0 / 2.0  )
                log_prob_0_1[2] = log(prob_zeta_mcmc[i, j])       + 2.0 * (par1 / 2.0 - 1.0) * log(bmat[j,j]) + (par1/ 2.0)*log(0.5) - loggamma(par1 / 2.0  )
                nzero_0_1[1] = sum_zero
                nzero_0_1[2] = sum_zero-1

                ### sparsa
                log_prob_0_1[1] += sum(like_acc[j:(k-1)])
                ### non sparsa
                z_prop[i,j] = 1.0
                from_qmat_to_qstar_triangularindex(i, j, qmat, qstar_prop, k, n, psi, z_prop)
                for jcol = j:(k-1)
                    compute_like_for_z(jcol, n, qstar_prop, like_prop, z_prop , k, ycentered)
                end
                for jcol = j:(k-1)
                    log_prob_0_1[2] += like_prop[jcol]
                end
                #log_prob_0_1[2] += sum(like_prop[j:(k-1)])

                #@assert !isnan(log_prob_0_1[2])
                if isnan(log_prob_0_1[2])
                    samp = 1
                else
                    samp = argmax(log_prob_0_1 .+ rand(Gumbel(0,1))) 
                end
                
                sum_zero = nzero_0_1[samp]
                z_mcmc[i,j] = Float64(samp-1)

                if z_mcmc[i,j] == 0.0 ## sparsa, vecchio valore
                    # ! colonna
                    #n_uno_vec_prop[i] = n_uno_vec[i]
                    #n_zero_vec_prop[i] = n_zero_vec[i]
                    # ! riga
                    n_uno_vec_prop[j] = n_uno_vec[j]
                    n_zero_vec_prop[j] = n_zero_vec[j]
                    
                    
                    z_prop[i,j] = 0.0
                    like_prop[j:(k-1)] = like_acc[j:(k-1)]
                    

                    qstar_prop[i,j] = qstar[i,j]
                    if (j+1) < i
                        for jcol  = (j+1):(i-1)
                            for irow = i:k
                                qstar_prop[irow,jcol] = qstar[irow,jcol]
                                #@assert irow<= k
                                #@assert jcol <= k
                                #@assert i>j
                            end
                        end
                    end
                    if i<k
                        for jcol = i:(k-1)                            
                            for irow = (jcol+1):k
                                qstar_prop[irow,jcol] = qstar[irow,jcol]
                            
                            end
                        end
                    end

                else 
                    # ! colonna
                    #n_uno_vec[i] = n_uno_vec_prop[i]
                    #n_zero_vec[i] = n_zero_vec_prop[i]
                    # ! riga
                    n_uno_vec[j] = n_uno_vec_prop[j]
                    n_zero_vec[j] = n_zero_vec_prop[j]

                    like_acc[j:(k-1)] = like_prop[j:(k-1)] 
                    #qstar .= qstar_prop
                    #for irow = (j+1):k
                    #    for jcol = j:(k-1)
                    #        qstar[irow,jcol] = qstar_prop[irow,jcol]
                    #    end
                    #end

                    qstar[i,j] = qstar_prop[i,j]
                    if (j+1) < i
                        for jcol  = (j+1):(i-1)
                            for irow = i:k
                                qstar[irow,jcol] = qstar_prop[irow,jcol]
                                
                            end
                        end
                    end
                    if i<k
                        for jcol = i:(k-1)                            
                            for irow = (jcol+1):k
                                qstar[irow,jcol] = qstar_prop[irow,jcol]
                                #@assert irow<= k
                                #@assert jcol <= k
                                #@assert i>j
                            end
                        end
                    end

                
                end


            else # zeta = 1

                # ! colonna
                #n_uno_vec_prop[i] = n_uno_vec[i] - 1.0
                #n_zero_vec_prop[i] = n_zero_vec[i] + 1.0
                #par0 = (nu - Float64(k) + 1.0  + maximum(n_uno_vec_prop) -n_uno_vec_prop[i] )
                #par1 = (nu - Float64(k) + 1.0  + maximum(n_uno_vec) -n_uno_vec[i] )
                # ! riga
                n_uno_vec_prop[j] = n_uno_vec[j] - 1.0
                n_zero_vec_prop[j] = n_zero_vec[j] + 1.0
                par0 = (nu - Float64(k) + 1.0  + n_uno_vec_prop[j] )
                par1 = (nu - Float64(k) + 1.0  + n_uno_vec[j] )
                
                

                log_prob_0_1[1] = log(1.0 - prob_zeta_mcmc[i, j]) + 2.0 * (par0 / 2.0 - 1.0) * log(bmat[j,j]) + (par0/ 2.0)*log(0.5) - loggamma(par0 / 2.0  )
                log_prob_0_1[2] = log(prob_zeta_mcmc[i, j])       + 2.0 * (par1 / 2.0 - 1.0) * log(bmat[j,j]) + (par1/ 2.0)*log(0.5) - loggamma(par1 / 2.0  )


                #log_prob_0_1[1] = log(1.0 - prob_zeta_mcmc[i, j]) + 2.0 * ((nu - maximum(n_zero_vec_prop) - n_uno_vec_prop[i]) / 2.0 - 1.0) * log(bmat[j, j]) + ((nu - maximum(n_zero_vec_prop) - n_uno_vec_prop[i]) / 2.0)*log(0.5) - loggamma( (nu - maximum(n_zero_vec_prop)- n_uno_vec_prop[i]) / 2.0 )
                #log_prob_0_1[2] = log(prob_zeta_mcmc[i, j]) + 2.0 * ((nu - maximum(n_zero_vec) - n_uno_vec[i]) / 2.0 - 1.0) * log(bmat[j, j]) + ((nu - maximum(n_zero_vec) - n_uno_vec[i]) / 2.0)*log(0.5) - loggamma( (nu - maximum(n_zero_vec)- n_uno_vec[i]) / 2.0 )
                
                nzero_0_1[2] = sum_zero
                nzero_0_1[1] = sum_zero+1


                ### non sparsa
                log_prob_0_1[2] += sum(like_acc[j:(k-1)])
                ### sparsa
                z_prop[i,j] = 0.0
                from_qmat_to_qstar_triangularindex(i, j, qmat, qstar_prop, k, n, psi, z_prop)
                for jcol = j:(k-1)
                    compute_like_for_z(jcol, n, qstar_prop, like_prop, z_prop,k, ycentered )
                end
                for jcol = j:(k-1)
                    log_prob_0_1[1] += like_prop[jcol]
                end
                #log_prob_0_1[1] += sum(like_prop[j:(k-1)])


                #@assert !isnan(log_prob_0_1[1])
                if isnan(log_prob_0_1[1])
                    samp = 2
                else
                    samp = argmax(log_prob_0_1 .+ rand(Gumbel(0,1))) 
                end
                
                sum_zero = nzero_0_1[samp]
                z_mcmc[i,j] = Float64(samp-1)

                if z_mcmc[i,j] == 0.0 ## sparsa, nuovo
                    
                    # ! colonna
                    #n_uno_vec[i] = n_uno_vec_prop[i]
                    #n_zero_vec[i] = n_zero_vec_prop[i]
                    # ! riga
                    n_uno_vec[j] = n_uno_vec_prop[j]
                    n_zero_vec[j] = n_zero_vec_prop[j]

                    like_acc[j:(k-1)] = like_prop[j:(k-1)] 
                    #qstar .= qstar_prop
                    #for irow = (j+1):k
                    #    for jcol = j:(k-1)
                    #        qstar[irow,jcol] = qstar_prop[irow,jcol]
                    #    end
                    #end


                    qstar[i,j] = qstar_prop[i,j]
                    if (j+1) < i
                        for jcol  = (j+1):(i-1)
                            for irow = i:k
                                qstar[irow,jcol] = qstar_prop[irow,jcol]
                                #@assert irow<= k
                                #@assert jcol <= k
                                #@assert i>j
                            end
                        end
                    end
                    if i<k
                        for jcol = i:(k-1)                            
                            for irow = (jcol+1):k
                                qstar[irow,jcol] = qstar_prop[irow,jcol]
                                #@assert irow<= k
                                #@assert jcol <= k
                                #@assert i>j
                            end
                        end
                    end

                else # noo sparsa, vecchio valore
                    # ! colonna
                    #n_uno_vec_prop[i] = n_uno_vec[i]
                    #n_zero_vec_prop[i] = n_zero_vec[i]
                    # ! riga
                    n_uno_vec_prop[j] = n_uno_vec[j]
                    n_zero_vec_prop[j] = n_zero_vec[j]

                    z_prop[i,j] = 1.0
                    like_prop[j:(k-1)] = like_acc[j:(k-1)]

                    #qstar_prop .= qstar
                    #for irow = (j+1):k
                    #    for jcol = j:(k-1)
                    #        qstar_prop[irow,jcol] = qstar[irow,jcol]
                    #    end
                    #end

                    qstar_prop[i,j] = qstar[i,j]
                    if (j+1) < i
                        for jcol  = (j+1):(i-1)
                            for irow = i:k
                                qstar_prop[irow,jcol] = qstar[irow,jcol]
                                #@assert irow<= k
                                #@assert jcol <= k
                                #@assert i>j
                            end
                        end
                    end
                    if i<k
                        for jcol = i:(k-1)                            
                            for irow = (jcol+1):k
                                qstar_prop[irow,jcol] = qstar[irow,jcol]
                                #@assert irow<= k
                                #@assert jcol <= k
                                #@assert i>j
                            end
                        end
                    end
                    
                end
                
            end

            

        end

    end
    from_bstar_to_b_mat(bstar_mcmc, b_mat_mcmc)
    from_bmat_to_all_q(b_mat_mcmc, q_mat_mcmc, q_mat_sparse_mcmc, z_mcmc, psi)
    lambda_mcmc[:, :] = q_mat_sparse_mcmc * transpose(q_mat_sparse_mcmc)
end


# ! se mai la cambio devo tenr conto de numeri degli 1 nella distribuzione delle chi quadro
#function sample_zeta_efficient!(iterMCMC::Int64, ydata::Matrix{Float64}, bstar_mcmc::Vector{Float64}, z_mcmc::Matrix{Float64}, model_sparse::T_SPARSE, type_mcmc_lambda::T_HMC, psi::Matrix{Float64}, nu::Float64, r_hmc_mcmc::Vector{Float64}, q_mat_mcmc::Matrix{Float64}, q_mat_sparse_mcmc::Matrix{Float64}, lambda_mcmc::Matrix{Float64}, b_mat_mcmc::Matrix{Float64}, ycentered::Matrix{Float64}, prob_zeta_mcmc::Matrix{Float64},  prec_q::Vector{Matrix{Float64}}) where {T_SPARSE<:TypeSparse,T_HMC<:TypeMCMC}

#    sum_zero::Int64 = abs(sum(z_mcmc .- 1)) + 1


#    #n_sample::Int64 = size(model_sparse.prob_vec_zeros)[1]
#    k::Int64 = size(z_mcmc, 1)
#    n::Int64 = size(ycentered, 2)
#    nm::Int64 = size(bstar_mcmc, 1) - k
#    hc::Int64 = 0
    

#    samp::Int64 = 0
#    z_prop::Matrix{Float64} = deepcopy(z_mcmc)

#    bmat::Matrix{Float64} = zeros(Float64, k, k)
#    # b contiene i parametri non constraint
#    from_bstar_to_bmat(bstar_mcmc, bmat, k, n)
#    qmat::Matrix{Float64} = zeros(Float64, k, k)
#    # qmat contiene psi*B
#    from_bmat_to_qmat(bmat, qmat, k, n, psi)
#    #qstar::Matrix{Float64} = zeros(Float64, k, k)
#    ## qmat è la matrice contrait
#    #from_qmat_to_qstar_init(qmat, qstar, k, n, psi, z_mcmc)


#    log_prob_0_1::Vector{Float64} = zeros(Float64, 2)
#    nzero_0_1::Vector{Int64} = zeros(Int64, 2)


#    lambda_algo::Matrix{Float64} = deepcopy(lambda_mcmc)
#    #lambda_algo_sim::Symmetric{Float64, Matrix{Float64}} = deepcopy(lambda_mcmc)
#    #lambda_algo::Symmetric{Float64, Matrix{Float64}} = deepcopy(lambda_mcmc)
#    lambda_nonconstr::Matrix{Float64} = qmat * transpose(qmat)

#    q_algo::Matrix{Float64} = zeros(Float64, k,k)
#    #q_nonconstr::Matrix{Float64} = zeros(Float64, k, k)
#    #b_nonconstr::Matrix{Float64} = zeros(Float64, k, k)
#    psi_inv::Matrix{Float64} = inv(psi)
#    psi_inv_t::Matrix{Float64} = transpose(psi_inv)
#    i::Int64 = k
#    j::Int64 = k-1

#    indexj::Int64 = 0
#    indexi::Int64 = 0
#    n_zero_vec = zeros(Float64, k)
#    for j = 1:(k-1)
#        for i = (j+1):k
#            n_zero_vec[j] += Float64((1.0 - z_mcmc[i, j]))
#        end
#    end
#    max_perm = Int64(trunc(k * (k - 1) / 2))
#    perm_vect::Vector{Int64} = zeros(Int64,k)
#    inv_perm_vect::Vector{Int64} = zeros(Int64, k)
#    for index_perm = 1:max_perm
        
#        perm_vect .= randperm(k)
#        inv_perm_vect .= sortperm(perm_vect)
#        if perm_vect[ k-1] > perm_vect[ k]
#            perm_vect[k-1], perm_vect[k] = perm_vect[k], perm_vect[k-1]
#        end
#        indexj = perm_vect[ k-1]
#        indexi = perm_vect[   k]
        

#        lambda_algo .= lambda_algo[:, perm_vect]
#        lambda_algo .= lambda_algo[perm_vect, :]
#        if isposdef(lambda_algo) == false
#            #println(isposdef(lambda_mcmc))
#            #println("a=", lambda_algo)
#            #println("b=", lambda_mcmc)
#            #println(index_perm)
#            #println(perm_mat[index_perm, :])
#            lambda_algo .= lambda_algo[:, inv_perm_vect]
#            lambda_algo .= lambda_algo[inv_perm_vect, :]
#            println("For break - iter ", iterMCMC)
#            break

#        end
#        q_algo .= cholesky(Symmetric(lambda_algo)).L

#        par_sd = q_algo[j, j]
#        par_reg_constr = 0.0
#        for ll in 1:(j-1)
#            par_reg_constr += q_algo[i, ll] * q_algo[j, ll]
#        end
#        par_reg_constr = -par_reg_constr / par_sd

#        par_reg_unconstr::Float64 = 0.0
#        m::Float64 = par_sd * psi[indexi, indexj] / psi[indexi, indexi]^2
#        v::Float64 = 1.0/psi[indexi, indexi]^2
#        if z_mcmc[indexi, indexj] == 0

#        elseif z_mcmc[indexi, indexj] == 1
#            for iobs = 1:n
#                v += ycentered[indexi, iobs]^2
#                m += -par_sd * ycentered[indexi, iobs]*ycentered[indexj, iobs]
#            end
#        end
#        v = 1.0 / v
#        m = v * m
#        par_reg_unconstr = rand(Normal(m, v^0.5))
        
#        ### z = 0
#        log_prob_0_1[1] = log(1.0 - prob_zeta_mcmc[i, j]) + 2.0 * ((nu - (0 + 1) - j + 1) / 2.0 - 1.0) * log(par_sd) + ((nu - (0 + 1) - j + 1) / 2.0) * log(0.5) - loggamma((nu - (0 + 1) - j + 1) / 2.0)
#        ### z = 1
#        log_prob_0_1[2] = log(prob_zeta_mcmc[i, j]) + 2.0 * ((nu - (0) - j + 1) / 2.0 - 1.0) * log(par_sd) + ((nu - (0) - j + 1) / 2.0) * log(0.5) - loggamma((nu - (0) - j + 1) / 2.0)
    
    
    
#        for iobs = 1:n
#            log_prob_0_1[1] += -0.5 * (ycentered[indexj, iobs] * par_sd + ycentered[indexi, iobs] * par_reg_constr)^2
#            log_prob_0_1[2] += -0.5 * (ycentered[indexj, iobs] * par_sd + ycentered[indexi, iobs] * par_reg_unconstr)^2
#        end
        
#        if isnan(log_prob_0_1[2])
#            samp = 1
#        else
#            samp = argmax(log_prob_0_1 .+ rand(Gumbel(0, 1)))
#        end

#        z_mcmc[indexi, indexj] = Float64(samp - 1)
#        if z_mcmc[indexi, indexj] == 0

#            q_algo[i, j] = par_reg_constr

#        elseif z_mcmc[indexi, indexj] == 1

#            q_algo[i, j] = par_reg_unconstr

#        end

#        lambda_algo .= q_algo * transpose(q_algo)
        
#        ### rimetto in ordine
#        lambda_algo .= lambda_algo[:, inv_perm_vect]
#        lambda_algo .= lambda_algo[inv_perm_vect, :]
#        #if isposdef(lambda_algo) == false
#        #    println("")
#        #    println("")
#        #    println("c=", q_algo)
            

#        #end

#    end
#    lambda_mcmc .= lambda_algo
#    q_constr::Matrix{Float64} = cholesky(Symmetric(lambda_mcmc)).L
#    q_nonconstr::Matrix{Float64} = deepcopy(q_constr)
#    index_1::Vector{Int64} = zeros(Int64,k)
#    index_0::Vector{Int64} = zeros(Int64, k)
#    num_zero::Int64 = 0
#    num_one::Int64 = 0
#    cond_var::Matrix{Float64} = zeros(Float64,k,k)
#    cond_mean::Vector{Float64} = zeros(Float64, k)
#    marg_mean::Vector{Float64} = zeros(Float64, k)
#    b_nonconstr::Matrix{Float64} = zeros(Float64, k, k)
#    for j = 1:(k-1)
        
#        num_zero = 0
#        num_one = 0
#        for i = (j+1):k
#            if z_mcmc[i, j] == 0
#                num_zero += 1
#                index_0[num_zero] = i
#            else
#                num_one += 1
#                index_1[num_one] = i
#            end
#        end
#        if num_zero > 0
#            if j != (k-1)
#                marg_mean .= psi[:, j] * q_constr[j, j]

#                cond_var[1:num_zero, 1:num_zero] .= inv(prec_q[j+1][index_0[1:num_zero].-j, index_0[1:num_zero].-j])


#                cond_mean[1:num_zero] .= marg_mean[index_0[1:num_zero]] - cond_var[1:num_zero, 1:num_zero] * prec_q[j+1][index_0[1:num_zero].-j, index_1[1:num_one].-j] * (q_constr[index_1[1:num_one]] - marg_mean[index_1[1:num_one]])

#                q_nonconstr[index_0[1:num_zero]] = rand(MvNormal(cond_mean[1:num_zero], Symmetric(cond_var[1:num_zero, 1:num_zero])))
#            else

#                q_nonconstr[k, k-1] = rand(Normal(psi[k, k-1] * q_nonconstr[k-1, k-1], psi[k, k]))
                

#            end
            
#        end
#    end

#    b_nonconstr .= psi_inv * q_nonconstr
#    ## j = 1
#    h::Int64 = 1
#    bstar_mcmc[h] = log(b_nonconstr[1,1])
#    h += 1
#    for i = 2:k
#        bstar_mcmc[h] = b_nonconstr[i, 1]
#        h += 1
#    end
#    for j = 2:(k-1)
#        bstar_mcmc[h] = log(b_nonconstr[j, j])
#        h += 1
#        for i = (j+1):k
#            bstar_mcmc[h] = b_nonconstr[i, j]
#            h += 1
#        end
#    end
#    bstar_mcmc[k] = log(b_nonconstr[k, k])
    
#    from_bstar_to_b_mat(bstar_mcmc, b_mat_mcmc)
#    from_bmat_to_all_q(b_mat_mcmc,  q_mat_mcmc, q_mat_sparse_mcmc, z_mcmc,psi)
#    lambda_mcmc[:,:] = q_mat_sparse_mcmc*transpose(q_mat_sparse_mcmc)
#end


