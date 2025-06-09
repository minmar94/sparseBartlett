
module SparseBartlettPackage

    using BandedMatrices
    using Distributions, Random
    using LinearAlgebra, PDMats, StatsBase
    using ToggleableAsserts
    using ProgressMeter
    using ForwardDiff
    using SpecialFunctions
    using Enzyme
    

    include(joinpath("types.jl"))
    include(joinpath("hmc.jl"))
    include(joinpath("functions.jl"))    
    include(joinpath("mcmc.jl"))
    include(joinpath("mcmc_glm.jl"))
    include(joinpath("det_prop.jl"))
    include(joinpath("sample_q.jl"))
    include(joinpath("sample_glm_q.jl"))
    include(joinpath("sample_z.jl"))
    #include(joinpath("sample_glm_z.jl"))
    include(joinpath("sample_mu.jl"))
    #include(joinpath("hmc_functions.jl"))
    

    export
        SparseRJ,
        GeneralS,
        HMCRJ,
        BaseHMC,
        GenSparse,
        HMC_DualAv_NoUTurn,
        HMC_DualAv,
        mcmc_glm,
        compute_unique_permutations,
        #StandardMCMC,
        #HMC,
        #TypeSparse,
        #NoSparse,
        #DiagonalS,
        #FixSparse,
        #SpikeSlabSparse,
        mcmc
end # module
