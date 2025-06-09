
function log_dens_h(hvec::Vector{Float64}, add_const::Float64)

    return  -0.5 * sum(hvec.^2.0)+ add_const
end

function leapfrog_nocheck(theta::Vector{Float64}, r_hmc::Vector{Float64}, epsilon::Float64, d_b::Vector{Float64}, log_dens::Function)
    
    
    gradient!(Reverse, d_b, log_dens, theta)
    r_hmc[:] .= r_hmc .+ (epsilon ./ 2.0) .* d_b
    
    theta[:] .= theta[:] .+ epsilon .* r_hmc
    
    gradient!(Reverse, d_b, log_dens, theta)
    r_hmc[:] .= r_hmc .+ (epsilon ./ 2.0) .* d_b
    
    return nothing
   
end

function leapfrog(theta::Vector{Float64}, r_hmc::Vector{Float64}, epsilon::Float64, d_b::Vector{Float64}, log_dens::Function)::Bool
    
    
    gradient!(Reverse, d_b, log_dens, theta)
    r_hmc[:] .= r_hmc + (epsilon / 2.0) .* d_b
    theta[:] .= theta[:] + epsilon * r_hmc
    gradient!(Reverse, d_b, log_dens, theta)
    r_hmc[:] .= r_hmc + (epsilon / 2.0) .* d_b
    
    return isnan(log_dens(theta)) || isinf(log_dens(theta)) || isnan(sum(r_hmc.^2.0)) || isinf(sum(r_hmc.^2.0)) || isnan(sum(d_b)) || isinf(sum(d_b)) 
   
end

function retry_leapfrog_type1(theta::Vector{Float64}, r_hmc::Vector{Float64}, theta_start::Vector{Float64}, r_hmc_start::Vector{Float64}, epsilon::Float64, d_b::Vector{Float64}, log_dens::Function)::Bool

    theta[:] .= theta_start[:]
    r_hmc[:] .= r_hmc_start[:]

    return leapfrog(theta, r_hmc, epsilon, d_b, log_dens)

end

function find_reasonable_epsilon(theta_start::Vector{Float64}, log_dens::Function, d_b::Vector{Float64}, epsilon::Vector{Float64}, add_const::Float64)
    
    iter_gen::Int64 = 0
    r_hmc_start::Vector{Float64} = rand(Normal(0.0,1.0),size(theta_start,1))
    
    theta = deepcopy(theta_start)
    r_hmc = deepcopy(r_hmc_start)
    #epsilon[1] = 1.0

    ltheta_start::Float64 = log_dens(theta)
    #lr_start::Float64 = sum(r_hmc.^2.0) 
    lr_start::Float64 = log_dens_h(r_hmc, add_const)
    fail::Bool = false


    epsilon_start::Float64 = epsilon[1]
    

    while fail
        iter_gen += 1
        if iter_gen >1
            epsilon[1] = epsilon[1]*0.5
            println("New test epsilon -  starting value ",  epsilon[1])
        end
        fail = leapfrog(theta, r_hmc, epsilon[1], d_b, log_dens)
    
        #println("A")
        iterator::Int64 = 0
        while fail
            
            epsilon[1] = epsilon[1] * 0.5
            
            fail =  retry_leapfrog_type1(theta, r_hmc, theta_start, r_hmc_start, epsilon[1], d_b, log_dens)
            
            #println([trunc.(theta, digits=3) ])
            iterator += 1
            if iterator > 100; error("iterator - find_reasonable_epsilon") end
        
        end
        #println("EPX ", epsilon[1])
        #println("B")
        iterator = 0



        a::Float64 = 2.0 * Int64( (log_dens(theta)+ log_dens_h(r_hmc, add_const)  - (ltheta_start +lr_start)) > log(0.5)) - 1.0
        
        #while (exp(log_dens(theta) - 0.5 * sum(r_hmc.^2.0)) / exp(ltheta_start - 0.5 * lr_start))^a > 2.0^(-a)

        
        while (a*(log_dens(theta) + log_dens_h(r_hmc, add_const)   - (ltheta_start + lr_start)) > -a * log(2.0) ) || fail
            
            theta_start[:] = theta[:]
            r_hmc_start[:] = r_hmc[:]
            
            epsilon[1] = 2.0^a * epsilon[1]
            println(epsilon[1])
            println([log_dens(theta), - 0.5 * sum(r_hmc.^2.0), - ltheta_start ,-lr_start])
            fail = leapfrog(theta, r_hmc, epsilon[1], d_b, log_dens)
        
            #while fail
            #    epsilon[1] = epsilon[1] * 0.5
            #    fail =  retry_leapfrog_type1(theta, r_hmc, theta_start, r_hmc_start, epsilon[1], d_b, log_dens)
                
            #    iterator += 1
            #    if iterator > 100; error("iterator - find_reasonable_epsilon step 2") end

            #end
        end

        if iter_gen >100
            error("iter_gen - find_reasonable_epsilon step 2")
        end
    end
    

    theta[:] = theta_start[:]
    
    #println("A ", sum(d_b))
    return nothing
end

function restart_hmc_par(type_hmc::HMC_DualAv)
    type_hmc.mu[1] = log(10.0*type_hmc.epsilon[1])
    type_hmc.bar_epsilon[:] = type_hmc.bar_epsilon_init[:]
    type_hmc.bar_h[:] = type_hmc.bar_h_init[:]
    type_hmc.gamma[:] = type_hmc.gamma_init[:]
    type_hmc.t[:] = type_hmc.t_init[:]
    type_hmc.k[:] = type_hmc.k_init[:]
end

function restart_hmc_par(type_hmc::HMC_DualAv_NoUTurn)
    type_hmc.mu[1] = log(10.0*type_hmc.epsilon[1])
    type_hmc.bar_epsilon[:] = type_hmc.bar_epsilon_init[:]
    type_hmc.bar_h[:] = type_hmc.bar_h_init[:]
    type_hmc.gamma[:] = type_hmc.gamma_init[:]
    type_hmc.t[:] = type_hmc.t_init[:]
    type_hmc.k[:] = type_hmc.k_init[:]
end


function dual_av_update_parameters(iter::Int64, alpha::Float64, type_hmc::T_HMC) where{T_HMC <: TypeMCMC}

    type_hmc.bar_h[1] = (1.0 - 1.0/(iter + type_hmc.t[1])) * type_hmc.bar_h[1] + 1.0/(iter + type_hmc.t[1])*(type_hmc.delta[1]-alpha)
    
    type_hmc.epsilon[1] = exp(type_hmc.mu[1] - sqrt(iter)/type_hmc.gamma[1]*type_hmc.bar_h[1])

    type_hmc.bar_epsilon[1] = exp( iter^(-type_hmc.k[1])*log(type_hmc.epsilon[1]) + (1.0- iter^(-type_hmc.k[1])) * log(type_hmc.bar_epsilon[1]) )

end

#### #### #### #### #### 
#### HAMILTONIANS
#### #### #### #### #### 
function hmc_sample(log_dens::Function, d_b::Vector{Float64},   theta::Vector{Float64},  r_hmc::Vector{Float64}, type_hmc::BaseHMC, iter::Int64, add_const::Float64) 
    base_hmc(log_dens, d_b,   theta,  r_hmc, type_hmc, iter, add_const)
end

#### ####  BASE HAMILTONIAN #### #### 
function base_hmc(log_dens::Function, d_b::Vector{Float64},   theta::Vector{Float64},  r_hmc::Vector{Float64}, type_hmc::BaseHMC, iter::Int64, add_const::Float64)

    #find_reasonable_epsilon(theta, log_dens, d_b, type_hmc.epsilon)


    npar::Int64 = size(theta,1)


    theta_start::Vector{Float64} = deepcopy(theta)
    r_hmc_start::Vector{Float64} = deepcopy(r_hmc)

    theta0::Vector{Float64} = deepcopy(theta_start)
    r_hmc0::Vector{Float64} = deepcopy(r_hmc)
    

    for m = 1:type_hmc.M[1]
        
        base_hmc_iter(m, log_dens, d_b,   theta,  r_hmc, type_hmc, iter, theta_start, r_hmc_start, theta0, r_hmc0, npar, add_const)

    end

    if type_hmc.restart[1]
        type_hmc.epsilon[1] = type_hmc.epsilon_init[1]
    end

end


function base_hmc_iter(m::Int64, log_dens::Function, d_b::Vector{Float64},   theta::Vector{Float64},  r_hmc::Vector{Float64}, type_hmc::T_HMC, iter::Int64, theta_start::Vector{Float64}, r_hmc_start::Vector{Float64}, theta0::Vector{Float64}, r_hmc0::Vector{Float64}, npar::Int64, add_const::Float64)::Float64 where{T_HMC <: TypeMCMC}
    
    iterator::Int64 = 0
    fail::Bool = false
    alpha_mh::Float64 = 0.0
    log_like_0::Float64 = 0.0

    r_hmc[:] = rand(Normal(0.0,1.0),npar)

    theta0[:] = theta[:]
    r_hmc0[:] = r_hmc[:]

    log_like_0 = log_dens(theta0) + log_dens_h(r_hmc0, add_const)  

    for i = 1:type_hmc.L[1]
        
        theta_start[:] .= theta[:]
        r_hmc_start[:] .= r_hmc[:]
        fail = leapfrog(theta, r_hmc, type_hmc.epsilon[1], d_b, log_dens)

        #iterator = 0
        #while fail

        #    type_hmc.epsilon[1] = type_hmc.epsilon[1] * 0.5
        #    fail =  retry_leapfrog_type1(theta, r_hmc, theta_start, r_hmc_start, type_hmc.epsilon[1], d_b, log_dens)
        #    iterator += 1

        #    if iterator >100; error("too many iterations - base_hmc") end

        #end

    end

    #####
    #alpha_mh = min(1.0, exp(log_like_0 - (log_dens(theta) - 0.5 * sum(r_hmc.^2.0))))
    if fail
        alpha_mh = 0.0
    else
        alpha_mh = min(1.0, exp(log_dens(theta) + log_dens_h(r_hmc, add_const)  - log_like_0 ) )
    end
    

    if rand(Uniform(0.0,1.0))< alpha_mh

        

    else
        
        theta[:] .= theta0[:]
        r_hmc[:] .= r_hmc0[:]

    end

    return alpha_mh

end


#### ####  HAMILTONIAN WITH DUAL AVERAGING #### #### 
function hmc_sample(log_dens::Function, d_b::Vector{Float64},   theta::Vector{Float64},  r_hmc::Vector{Float64}, type_hmc::HMC_DualAv, iter::Int64, add_const::Float64) 
    hmc_dual_av(log_dens, d_b,   theta,  r_hmc, type_hmc, iter, add_const)
end


function hmc_dual_av(log_dens::Function, d_b::Vector{Float64},   theta::Vector{Float64},  r_hmc::Vector{Float64}, type_hmc::HMC_DualAv, iter::Int64, add_const::Float64)

    #
    if type_hmc.restart[1]
        restart_hmc_par(type_hmc)
    end
    if type_hmc.general_iter[1] == 1
        
        find_reasonable_epsilon(theta, log_dens, d_b, type_hmc.epsilon,  add_const)
        type_hmc.mu[1] = log(10.0*type_hmc.epsilon[1])
        
    end

    npar::Int64 = size(theta,1)

    theta_start::Vector{Float64} = deepcopy(theta)
    r_hmc_start::Vector{Float64} = deepcopy(r_hmc)

    theta0::Vector{Float64} = deepcopy(theta_start)
    r_hmc0::Vector{Float64} = deepcopy(r_hmc)
    alpha::Float64 = 0.0
    for m = 1:type_hmc.M[1]
        
        type_hmc.general_iter[1] += 1
        type_hmc.L[1] = max(1, round(type_hmc.lambda[1]/type_hmc.epsilon[1]))
        type_hmc.L[1] = min(100, type_hmc.L[1])

        
        alpha = base_hmc_iter(m, log_dens, d_b,   theta,  r_hmc, type_hmc, iter, theta_start, r_hmc_start, theta0, r_hmc0, npar, add_const)

        if iter < type_hmc.iter_adapt[1]
            dual_av_update_parameters(type_hmc.general_iter[1], alpha, type_hmc)
        end
    end

end



#### ####  HAMILTONIAN WITH DUAL AVERAGING and NOUTURN #### #### 
function build_tree_hmc(log_dens::Function, d_b::Vector{Float64},  utree::Float64, vtree::Int64, jtree::Int64, epsilon::Float64, deltamax::Float64, theta::Vector{Float64},  r_hmc::Vector{Float64}, theta_zero::Vector{Float64},  r_hmc_zero::Vector{Float64}, add_const::Float64)

    
    
    r_hmc_prime = deepcopy(r_hmc)
    theta_prime = deepcopy(theta)

    total_dens::Float64 = 0.0
    #ltheta_start::Float64 = log_dens(theta)
    #lr_start::Float64 = sum(r_hmc.^2.0)

    ntree_prime::Int64 = 0
    stree_prime::Int64 = 0
    ntree_prime_v2::Int64 = 0
    stree_prime_v2::Int64 = 0
    #println("eps ", epsilon, " - ", jtree)
    if jtree == 0
        # Base case - take one leapfrog step in the direction v.
        leapfrog(theta_prime, r_hmc_prime, epsilon, d_b, log_dens)
        # NOTE: this trick prevents the log-joint or its graident from being infinte
        
        iterator::Int64 = 0
        while isnan(log_dens(theta_prime)) || isinf(log_dens(theta_prime)) || isnan(sum(r_hmc_prime.^2.0)) || isinf(sum(r_hmc_prime.^2.0))

            epsilon = epsilon * 0.5
            theta_prime[:] .= theta[:]
            r_hmc_prime[:] .= r_hmc[:]
            
            leapfrog(theta_prime, r_hmc_prime, epsilon, d_b, log_dens)
            iterator += 1

            if iterator >100
                error("iterator v2")
            end

        end
        #println("jtree=0 ", epsilon)
        total_dens = log_dens(theta_prime) + log_dens_h(r_hmc_prime, add_const) 
        ntree_prime = Int64(utree <= exp(total_dens))
        stree_prime = Int64(utree < exp(deltamax + total_dens))

        #println("CHECK")
        #println([total_dens , log_dens(theta_zero), - 0.5 * sum(r_hmc_zero.^2.0) ,  log_dens(theta_prime), - 0.5 * sum(r_hmc_prime.^2.0)])
        return theta_prime, r_hmc_prime, theta_prime, r_hmc_prime, theta_prime, ntree_prime, stree_prime, min(1.0, exp(total_dens - (log_dens(theta_zero) + log_dens_h(r_hmc_zero, add_const)   ))), 1.0

    else

        theta_minus, r_hmc_minus, theta_plus, r_hmc_plus,  theta_prime, ntree_prime, stree_prime, a_prime, nalpha_prime = build_tree_hmc(log_dens, d_b,  utree, vtree, jtree-1, epsilon, deltamax, theta,  r_hmc, theta_zero,  r_hmc_zero, add_const)

        if stree_prime == 1
            if vtree == -1

                theta_minus, r_hmc_minus,_,_,  theta_prime_v2, ntree_prime_v2, stree_prime_v2, a_prime_v2, nalpha_prime_v2 = build_tree_hmc(log_dens, d_b,  utree, vtree, jtree-1, epsilon, deltamax, theta_minus,  r_hmc_minus, theta_zero,  r_hmc_zero, add_const)
        
                

            else

                _, _, theta_plus, r_hmc_plus,  theta_prime_v2, ntree_prime_v2, stree_prime_v2, a_prime_v2, nalpha_prime_v2 = build_tree_hmc(log_dens, d_b,  utree, vtree, jtree-1, epsilon, deltamax, theta_plus,  r_hmc_plus, theta_zero,  r_hmc_zero, add_const)

                

            end
            if rand(Uniform(0.0,1.0)) < ntree_prime_v2/(ntree_prime_v2 + ntree_prime)
                theta_prime[:] = theta_prime_v2[:]
            end
            
            stree_prime = stree_prime_v2 * Int64( sum((theta_plus .- theta_minus).*r_hmc_minus)>=0.0 ) * Int64( sum((theta_plus .- theta_minus).*r_hmc_plus)>=0.0 )

            #println("B ", [stree_prime_v2, Int64( sum((theta_plus .- theta_minus).*r_hmc_minus)>=0.0 ), Int64( sum((theta_plus .- theta_minus).*r_hmc_plus)>=0.0 )])

            ntree_prime = ntree_prime_v2 + ntree_prime
            
            a_prime = a_prime + a_prime_v2
            nalpha_prime = nalpha_prime + nalpha_prime_v2

        end
        return theta_minus, r_hmc_minus, theta_plus, r_hmc_plus,  theta_prime, ntree_prime, stree_prime, a_prime, nalpha_prime
    end
end


function hmc_sample(log_dens::Function, d_b::Vector{Float64},   theta::Vector{Float64},  r_hmc::Vector{Float64}, type_hmc::HMC_DualAv_NoUTurn, iter::Int64, add_const::Float64) 
    hmc_dual_av_nouturn(log_dens, d_b,   theta,  r_hmc, type_hmc, iter, add_const)
end


function hmc_dual_av_nouturn(log_dens::Function, d_b::Vector{Float64},   theta::Vector{Float64},  r_hmc::Vector{Float64}, type_hmc::HMC_DualAv_NoUTurn, iter::Int64, add_const::Float64)

    #
    if type_hmc.restart[1]
        restart_hmc_par(type_hmc)
    end
    if type_hmc.general_iter[1] == 1

        find_reasonable_epsilon(theta, log_dens, d_b, type_hmc.epsilon, add_const)
        type_hmc.mu[1] = log(10.0*type_hmc.epsilon[1])
    end

    theta_start::Vector{Float64} = deepcopy(theta)
    r_hmc_start::Vector{Float64} = deepcopy(r_hmc)

    theta0::Vector{Float64} = deepcopy(theta_start)
    r_hmc0::Vector{Float64} = deepcopy(r_hmc)
    alpha::Float64 = 0.0

    theta_minus::Vector{Float64}  = deepcopy(theta)
    theta_plus::Vector{Float64}  = deepcopy(theta)
    

    r_hmc_minus::Vector{Float64}  = deepcopy(r_hmc)
    r_hmc_plus::Vector{Float64}  = deepcopy(r_hmc)
    


    for m = 1:type_hmc.M[1]

        type_hmc.general_iter[1] += 1
        alpha = hmc_dual_av_nouturn_iter(m, log_dens, d_b,   theta,  r_hmc, type_hmc, theta_minus, theta_plus, r_hmc_minus, r_hmc_plus, theta_start, r_hmc_start, add_const)
        
        if iter < type_hmc.iter_adapt[1]
            dual_av_update_parameters(type_hmc.general_iter[1], alpha, type_hmc)
        end
    end

end



function hmc_dual_av_nouturn_iter(m::Int64, log_dens::Function, d_b::Vector{Float64},   theta::Vector{Float64},  r_hmc::Vector{Float64}, type_hmc::HMC_DualAv_NoUTurn,  theta_minus::Vector{Float64}, theta_plus::Vector{Float64}, r_hmc_minus::Vector{Float64}, r_hmc_plus::Vector{Float64}, theta_start::Vector{Float64}, r_hmc_start::Vector{Float64}, add_const::Float64)

    
    r_hmc[:] = rand(Normal(0.0,1.0), size(theta,1))

    theta_zero = deepcopy(theta)
    r_hmc_zero = deepcopy(r_hmc)
    

    #add_const::Float64 = -min(log_dens(theta_zero) , log_dens_h(r_hmc_zero, 0.0) )
    add_const::Float64 = - log_dens(theta_zero) - log_dens_h(r_hmc_zero, 0.0) 

    
    utree = rand(Uniform(0.0, exp(log_dens(theta_zero) +log_dens_h(r_hmc_zero, add_const)  )))
    #println([utree, exp(log_dens(theta_zero) +log_dens_h(r_hmc_zero, add_const)  )])
    jtree::Int64 = 0
    ntree::Int64 = 1
    stree::Int64 = 1

    
    
    ntree_prime::Int64 = 1
    stree_prime::Int64 = 1
    vtree::Int64 = 1
    a::Float64 = 1
    nalpha::Int64 = 1

    theta_minus[:] = theta_zero[:]
    theta_plus[:] = theta_zero[:]
    theta_start[:] = theta_zero[:]

    r_hmc_minus[:] = r_hmc_zero[:]
    r_hmc_plus[:] = r_hmc_zero[:]
    r_hmc_start[:] = r_hmc_zero[:]


    while stree == 1
        
        vtree = rand([-1, 1])

        if vtree == -1

            theta_minus, r_hmc_minus,_,_,  theta_prime, ntree_prime, stree_prime, a, nalpha = build_tree_hmc(log_dens, d_b,  utree, vtree, jtree, type_hmc.epsilon[1], type_hmc.deltamax[1], theta_minus,  r_hmc_minus, theta_zero, r_hmc_zero, add_const)
    
            

        else

            _, _, theta_plus, r_hmc_plus,  theta_prime, ntree_prime, stree_prime, a, nalpha = build_tree_hmc(log_dens, d_b,  utree, vtree, jtree, type_hmc.epsilon[1], type_hmc.deltamax[1], theta_plus,  r_hmc_plus, theta_zero, r_hmc_zero, add_const)

            

        end
        
        if stree_prime == 1
            #println("stree_prime ", ntree_prime/ntree)
            if rand(Uniform(0.0,1.0)) < min(1.0,  ntree_prime/ntree)
                theta[:] = theta_prime[:]
            end

        end
        ntree = ntree + ntree_prime
        stree= stree_prime * Int64( sum((theta_plus .- theta_minus).*r_hmc_minus) >=0.0 ) * Int64( sum((theta_plus .- theta_minus).*r_hmc_plus)>=0.0 )
        jtree += 1
    end

    return Float64(a/nalpha)

end
