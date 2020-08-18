module StochasticProcesses

using Statistics
using Distributions

# Simulate pure birth process trajectory
function sim_traj_birth_ct(λ::Float64,tmax::Float64)::NTuple{2,Vector{Float64}}

    ## ME
    # ∂/∂t p(X(t)=x) = -λ⋅p(X(t)=x) + λ⋅p(X(t)=x-1)
    ## Solution:
    # X(t) ~ Poisson(λt)

    # Gillespie until reach tmax
    t = 0.0
    tvec = [t]
    while t<=tmax
        
        # Increase time counter
        t += rand(Exponential(1.0/λ))
        
        # Store time
        push!(tvec,t)

    end

    # Return (t,x(t)) tuple
    return (tvec,0:(length(tvec)-1))

end

# Simulate linear birth process trajectory
function sim_traj_birth_lin(λ::Float64,tmax::Float64)::NTuple{2,Vector{Float64}}

    ## ME: 
    # ∂/∂t p(X(t)=x) = -λ⋅x⋅p(X(t)=x) + λ⋅(x-1)⋅p(X(t)=x+1)
    ## Solution:
    # X(t) ~ TODO

    # Gillespie until reach X(t)=0
    x = 1.0
    t = 0.0
    tvec = [t]
    while t<=tmax
        
        # Increase time counter
        t += rand(Exponential(1.0/(λ*x)))
        
        # Store time
        push!(tvec,t)

        # Decrease alive counter
        x += 1.0

    end

    # Return (t,x(t)) tuple
    return (tvec,0:(length(tvec)-1))

end

# Simulate linear death process trajectory
function sim_traj_death_ct(μ::Float64,x0::Int64)::NTuple{2,Vector{Float64}}

    ## ME: 
    # ∂/∂t p(X(t)=x) = -μ⋅p(X(t)=x) + μ⋅p(X(t)=x+1)
    ## Solution (induction):
    # P[X(t)=x0+j|X0=x0] = x0⋅(x0+1)⋯(x0+j-1)⋅exp{-λ⋅x0⋅t}⋅(1-exp{-λ⋅t})^{j} / j

    # Gillespie until reach X(t)=0
    x = x0
    t = 0.0
    tvec = [t]
    while x>0
        
        # Increase time counter
        t += rand(Exponential(1.0/μ))
        
        # Store time
        push!(tvec,t)

        # Decrease alive counter
        x -= 1 

    end

    # Return (t,x(t)) tuple
    return (tvec,x0:-1:0)

end

# Simulate linear death process trajectory
function sim_traj_death_lin(μ::Float64,x0::Int64)::NTuple{2,Vector{Float64}}

    ## ME: 
    # ∂/∂t p(X(t)=x) = -μ⋅x⋅p(X(t)=x) + μ⋅(x+1)⋅p(X(t)=x+1)
    ## Solution:
    # X(t) ~ Binomial(exp(-μt);X0)

    # Gillespie until reach X(t)=0
    x = x0
    t = 0.0
    tvec = [t]
    while x>0
        
        # Increase time counter
        t += rand(Exponential(1.0/(μ*x)))
        
        # Store time
        push!(tvec,t)

        # Decrease alive counter
        x -= 1 

    end

    # Return (t,x(t)) tuple
    return (tvec,x0:-1:0)

end

# Simulate linear birth-death process trajectory
function sim_traj_birth_death_lin(λ::Float64,μ::Float64,tmax::Float64)::NTuple{2,Vector{Float64}}

    ## ME: 
    # ∂/∂t p(X(t)=x) = λ⋅(x-1)⋅p(X(t)=x-1) + μ⋅(x+1)⋅p(X(t)=x+1) - (λ+μ)⋅x⋅p(X(t)=x)
    ## Solution:
    # X(t) ~ 

    # Gillespie until reach X(t)=0
    xvec = [1]
    tvec = [0.0]
    while tvec[end]<tmax
        
        # Increase time counter
        if xvec[end]==0
            Δt = tmax
        else
            Δt = rand(Exponential(abs(1.0/((λ+μ)*xvec[end]))))
        end

        # Store time
        push!(tvec,tvec[end]+Δt)

        # Which reaction happens
        if xvec[end]==0
            Δx = 0
        elseif xvec[end]==1
            Δx = 1
        else
            Δx = rand()>=λ/(λ+μ) ? -1 : 1
        end

        # Decrease alive counter
        push!(xvec,xvec[end]+Δx)

    end

    # Return (t,x(t)) tuple
    return (tvec,xvec)

end

end # module
