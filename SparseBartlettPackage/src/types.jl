
abstract type TypeMCMC end   
struct StandardMCMC <: TypeMCMC

    function StandardMCMC()

        new()

    end
     
end

struct HMC <: TypeMCMC

    epsilon::Vector{Float64}
    iter::Vector{Int64}
    function HMC(epsilon::Float64, iter::Int64)

        new([epsilon], [iter])

    end
     
end

struct HMCRJ <: TypeMCMC

    epsilon::Vector{Float64}
    iter::Vector{Int64}
    delta::Vector{Float64}
    m_max::Vector{Int64}
    function HMCRJ(epsilon::Float64, iter::Int64, delta::Float64, m_max::Int64)

        new([epsilon], [iter], [delta], [m_max])

    end
     
end


struct BaseHMC <: TypeMCMC

    epsilon::Vector{Float64}
    L::Vector{Int64}
    M::Vector{Int64}
    restart::Vector{Bool}
    epsilon_init::Vector{Float64}
    L_init::Vector{Int64}
    M_init::Vector{Int64}
    general_iter::Vector{Int64}
    function BaseHMC(;epsilon::Float64, L::Int64, M::Int64, restart::Bool)

        new([epsilon], [L], [M], [restart], [epsilon], [L], [M], [1])

    end
     
end

struct HMC_DualAv <: TypeMCMC

    epsilon::Vector{Float64}
    L::Vector{Int64}
    M::Vector{Int64}
    delta::Vector{Float64}
    lambda::Vector{Float64}
    iter_adapt::Vector{Int64}

    restart::Vector{Bool}
    epsilon_init::Vector{Float64}
    L_init::Vector{Int64}
    M_init::Vector{Int64}
    general_iter::Vector{Int64}
    delta_init::Vector{Float64}
    lambda_init::Vector{Float64}

    mu::Vector{Float64}
    bar_epsilon::Vector{Float64}
    bar_h::Vector{Float64}
    gamma::Vector{Float64}
    t::Vector{Float64}
    k::Vector{Float64}

    mu_init::Vector{Float64}
    bar_epsilon_init::Vector{Float64}
    bar_h_init::Vector{Float64}
    gamma_init::Vector{Float64}
    t_init::Vector{Float64}
    k_init::Vector{Float64}

    function HMC_DualAv(;epsilon::Float64, L::Int64, M::Int64, restart::Bool, delta::Float64, lambda::Float64, iter_adapt::Int64   )

        new([epsilon], [L], [M], [delta], [lambda], [iter_adapt], [restart], [epsilon], [L], [M], [1], [delta], [lambda], [log(10.0)], [1.0], [0.0], [0.5], [10.0], [0.75], [log(10.0)], [1.0], [0.0], [0.5], [10.0], [0.75])

    end
     
end

struct HMC_DualAv_NoUTurn <: TypeMCMC

    epsilon::Vector{Float64}
    L::Vector{Int64}
    M::Vector{Int64}
    delta::Vector{Float64}
    #lambda::Vector{Float64}
    iter_adapt::Vector{Int64}

    restart::Vector{Bool}
    epsilon_init::Vector{Float64}
    L_init::Vector{Int64}
    M_init::Vector{Int64}
    general_iter::Vector{Int64}
    delta_init::Vector{Float64}
    #lambda_init::Vector{Float64}

    mu::Vector{Float64}
    bar_epsilon::Vector{Float64}
    bar_h::Vector{Float64}
    gamma::Vector{Float64}
    t::Vector{Float64}
    k::Vector{Float64}

    mu_init::Vector{Float64}
    bar_epsilon_init::Vector{Float64}
    bar_h_init::Vector{Float64}
    gamma_init::Vector{Float64}
    t_init::Vector{Float64}
    k_init::Vector{Float64}
    deltamax::Vector{Float64}

    function HMC_DualAv_NoUTurn(;epsilon::Float64, L::Int64, M::Int64, restart::Bool, delta::Float64,  iter_adapt::Int64   )

        new([epsilon], [L], [M], [delta], [iter_adapt], [restart], [epsilon], [L], [M], [1], [delta],  [log(10.0)], [1.0], [0.0], [0.5], [10.0], [0.75], [log(10.0)], [1.0], [0.0], [0.5], [10.0], [0.75], [1000.0])

    end
     
end
#########


abstract type TypeSparse end   
struct NoSparse <: TypeSparse

    function NoSparse()

        new()

    end
     
end
struct GenSparse <: TypeSparse

    prob_vec_zeros::Vector{Float64}
    function GenSparse(prob_vec::Vector{Float64})

        new(deepcopy(prob_vec))

    end
     
end

struct FixSparse <: TypeSparse

    zero_indices::Matrix{Int64}

    function FixSparse(zero_indices::Matrix{Int64})

        new(zero_indices)

    end
     
end

struct SpikeSlabSparse <: TypeSparse

    tau_small::Float64
    tau_large::Float64
    z_prior::Binomial{Float64}

    function SpikeSlabSparse(tau_small::Float64, tau_large::Float64, z_prior::Binomial{Float64})

        new(tau_small, tau_large, z_prior)

    end
     
end
struct SparseRJ <: TypeSparse

    prob::Vector{Float64}

    function SparseRJ(prob::Vector{Float64})

        new(prob)

    end
     
end

abstract type TypeS end   
struct DiagonalS <: TypeS

    function DiagonalS()

        new()

    end
     
end


struct GeneralS <: TypeS

    function GeneralS()

        new()

    end
     
end


#############



