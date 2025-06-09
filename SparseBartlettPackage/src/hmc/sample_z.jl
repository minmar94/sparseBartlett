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
    for s = 1:nm

        i = sample(1:k,1)[1]
        j = sample(1:k,1)[1]

        if i != j

            j,i = min(i,j), max(i,j)
            if z_mcmc[i,j] == 0.0
    
                #log_prob_0_1[1] = model_sparse.prob_vec_zeros[sum_zero] + log(1.0-prob_zeta_mcmc[i,j])
                #log_prob_0_1[2] = model_sparse.prob_vec_zeros[sum_zero-1]  + log(prob_zeta_mcmc[i,j])
                log_prob_0_1[1] = log(1.0-prob_zeta_mcmc[i,j])
                log_prob_0_1[2] = log(prob_zeta_mcmc[i,j])
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
                    
                    z_prop[i,j] = 0.0
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
                            
                            end
                        end
                    end

                else 
                    
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


            else
                #println([log(prob_zeta_mcmc[i,j]) , log(sum_prob_over1 + prob_zeta_mcmc[i,j] )])
                #log_prob_0_1[2] = model_sparse.prob_vec_zeros[sum_zero]+ log(prob_zeta_mcmc[i,j])
                #log_prob_0_1[1] = model_sparse.prob_vec_zeros[sum_zero+1] + log(1.0-prob_zeta_mcmc[i,j])
                log_prob_0_1[2] = log(prob_zeta_mcmc[i,j])
                log_prob_0_1[1] = log(1.0-prob_zeta_mcmc[i,j])
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

end




#function sample_zeta_and_q!(iterMCMC::Int64, ydata::Matrix{Float64}, bstar_mcmc::Vector{Float64}, z_mcmc::Matrix{Float64}, model_sparse::T_SPARSE,  type_mcmc_lambda::T_HMC, psi::Matrix{Float64}, nu::Float64, r_hmc_mcmc::Vector{Float64}, q_mat_mcmc::Matrix{Float64}, q_mat_sparse_mcmc::Matrix{Float64}, lambda_mcmc::Matrix{Float64}, b_mat_mcmc::Matrix{Float64}, psi_inv::Matrix{Float64}, array_par_inv::Array{Float64,3})  where {T_SPARSE<:TypeSparse, T_HMC<:TypeMCMC}

#    k::Int64 = size(z_mcmc,1)
#    nm::Int64 = size(bstar_mcmc,1)- k
#    n::Int64 = size(ydata,2)

#    qmat::Matrix{Float64} = zeros(Float64, k,k) 
#    bmat::Matrix{Float64} = zeros(Float64, k,k) 
#    qstar::Matrix{Float64} = zeros(Float64, k,k) 
#    ###### Creo la matrice B
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
#            qmat[i,j] = sum(psi[i,j:k] .* bmat[j:k,j])
#        end
#    end

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

#    z_orig::Int64 = 0
#    q_orig::Float64 = 0.0

#    sum_zero::Int64 = abs(sum(z_mcmc  .- 1)) +  1
    


#    i::Int64 = 1
#    j::Int64 = 1

#    mu_marg::Float64 = 0.0
#    mu_const::Float64 = 0.0
#    q_prop::Float64 = 0.0
#    z_prop::Int64 = 0
    

#    log_prob_prop::Float64 = 0.0
#    log_prob_acc::Float64 = 0.0

#    MH::Float64 = 0.0
    
#    vec_molt_prec::Vector{Float64} = zeros(Float64,k)
#    for s = 1:(2*nm)

#        i = sample(1:k,1)[1]
#        j = sample(1:k,1)[1]

#        if i != j

#            j,i = min(i,j), max(i,j)


#            z_orig = z_mcmc[i,j]
#            q_orig = qmat[i,j]
#            vec_molt_prec[:] = -array_par_inv[j,i,:] ./ array_par_inv[j,i,i]
#            vec_molt_prec[i] = 0.0

            
#            mu_marg = bmat[j,j]*psi[i,j] + sum(vec_molt_prec .* (qmat[:,j] .- bmat[j,:]))
        
#            mu_const = - sum(qstar[i,1:(j-1)] .* qstar[j,1:(j-1)]) ./ qstar[j,j]

#            q_prop = rand(Normal( (1-z_mcmc[i,j]) * mu_const + z_mcmc[i,j] * mu_marg, 1.0 ))
#            log_prob_prop = logpdf(Normal( (1-z_mcmc[i,j]) * mu_const + z_mcmc[i,j] * mu_marg, 1.0 ), q_prop)

#            z_prop = sample(0:1,1)[1] 
#            log_prob_acc = logpdf(Normal( (1-z_prop) * mu_const + z_prop * mu_marg, 1.0 ), qmat[i,j])

#            ### MH proposal
#            MH += log_prob_acc - log_prob_prop
#            ### MH prior
#            MH += logpdf(Normal(mu_marg, (1.0/array_par_inv[j,i,i])^0.5), q_prop)
#            MH -= logpdf(Normal(mu_marg, (1.0/array_par_inv[j,i,i])^0.5), qmat[i,j])
            
                
#            MH -= model_sparse.prob_vec_zeros[sum_zero]
#            if z_orig == 0
#                if z_prop == 0
#                    MH += model_sparse.prob_vec_zeros[sum_zero] 
#                else
#                    MH += model_sparse.prob_vec_zeros[sum_zero-1]
#                end
#            end
#            if z_orig == 1
#                if z_prop == 0
#                    MH += model_sparse.prob_vec_zeros[sum_zero+1] 
#                else
#                    MH += model_sparse.prob_vec_zeros[sum_zero]
#                end
#            end

#            #### Like accepted -  solo i futuri
#            for js = j:(k-1)
#                for is = (js+1):k
#                    qstar[is,js] = z_mcmc[is,js] * qmat[is,js] + (1.0-z_mcmc[is,js]) * (- sum(qstar[is,1:(js-1)] .* qstar[js,1:(js-1)]) / qstar[js,js])
#                end
#            end
#            for iobs = 1:n
#                for js = j:k
#                    MH -=  -0.5*( sum(ydata[js:k,iobs] .* qstar[js:k,js] ))^2.0
#                end
#            end

#            #### Like prop -  solo i futuri
#            z_mcmc[i,j] = z_prop
#            qmat[i,j] = q_prop
#            for js = j:(k-1)
#                for is = (js+1):k
#                    qstar[is,js] = z_mcmc[is,js] * qmat[is,js] + (1.0-z_mcmc[is,js]) * (- sum(qstar[is,1:(js-1)] .* qstar[js,1:(js-1)]) / qstar[js,js])
#                end
#            end
#            for iobs = 1:n
#                for js = j:k
#                    MH +=  -0.5*( sum(ydata[js:k,iobs] .* qstar[js:k,js] ))^2.0
#                end
#            end
            
#            if rand(Uniform(0.0,1.0))< exp(MH)
#                if z_orig != z_prop
#                    println("From ", z_orig, " to ", z_prop)
#                end
#                if z_orig == 0
#                    if z_prop == 1
#                        sum_zero +=  -1
#                    end
#                end
#                if z_orig == 1
#                    if z_prop == 0
#                        sum_zero +=  1
#                    end
#                end

#            else
#                z_mcmc[i,j] = z_orig
#                qmat[i,j] = q_orig
#                for js = j:(k-1)
#                    for is = (js+1):k
#                        qstar[is,js] = z_mcmc[is,js] * qmat[is,js] + (1.0-z_mcmc[is,js]) * (- sum(qstar[is,1:(js-1)] .* qstar[js,1:(js-1)]) / qstar[js,js])
#                    end
#                end
#            end


#        end

#    end
#    bmat[:,:] = psi_inv*qmat

#    h = 1
#    bstar_mcmc[h] = log(bmat[1,1])
#    h += 1
#    for i = 2:k
#        bstar_mcmc[h] = bmat[i,1]
#        h += 1
#    end
#    for j = 2:(k-1)
#        bstar_mcmc[h] = log(bmat[j,j])
#        h += 1
#        for i = (j+1):k
#            bstar_mcmc[h] = bmat[i,j]
#            h += 1
#        end
#    end
#    bstar_mcmc[h] = log(bmat[k,k])

#end