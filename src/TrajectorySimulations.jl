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

# Simulate birth-death process trajectory (closed system)
function sim_traj_birth_death_lin_closed(κ1::Float64,κ2::Float64,ntot::Int64,tmax::Float64)::NTuple{2,Vector{Float64}}

    ## Reaction network: 
    # X1 → X2 (a particle becomes alive)
    # X1 ← X2 (a particle becomes dead)
    ## ME:
    # ∂/∂t p(X(t)=[x1,x2]) = κ1 ⋅ (x1+1) ⋅p(X(t)=[x1+1,x2-1]) 
    #                       + κ2 ⋅ (x2+1) ⋅ p(X(t)=[x1-1,x2+1]) 
    #                       - (κ1⋅x1+κ2⋅x2) ⋅ p(X(t)=[x1,x2])

    # Gillespie
    tvec = [0.0]
    x2vec = [ntot]
    while tvec[end]<tmax
        
        # Sum propensity functions
        a = κ1*(ntot-x2vec[end])+κ2*x2vec[end]

        # Time increment
        Δt = rand(Exponential(1.0/a))

        # Store time
        push!(tvec,tvec[end]+Δt)

        # Which reaction happens
        Δx = rand()>=(κ1*(ntot-x2vec[end])/a) ? -1 : +1

        # Decrease alive counter
        push!(x2vec,x2vec[end]+Δx)

    end

    # Return (t,x(t)) tuple
    return (tvec,x2vec)

end

# Simulate birth-death process trajectory (open system)
function sim_traj_birth_death_lin_open(κ1::Float64,κ2::Float64,tmax::Float64)::NTuple{2,Vector{Float64}}

    ## Reaction network: 
    # ∅ → X2 (a particle becomes alive)
    # X1 ← ∅ (a particle becomes dead)
    ## ME:
    # ∂/∂t p(X(t)=x) = κ1 ⋅ p(X(t)=x-1) + κ2 ⋅ (x+1) ⋅ p(X(t)=x+1) - (κ1+κ2⋅x) ⋅ p(X(t)=x)

    # Gillespie
    xvec = [0]
    tvec = [0.0]
    while tvec[end]<tmax
        
        # Sum propensity functions
        a = κ1+κ2*xvec[end]

        # Time increment
        Δt = rand(Exponential(1.0/a))

        # Store time
        push!(tvec,tvec[end]+Δt)

        # Which reaction happens
        Δx = rand()>=(κ1/a) ? -1 : +1

        # Decrease alive counter
        push!(xvec,xvec[end]+Δx)

    end

    # Return (t,x(t)) tuple
    return (tvec,xvec)

end

# Simulate SIR model
function sim_traj_sir(κ1::Float64,κ2::Float64,ntot::Int64,tmax::Float64)::NTuple{4,Vector{Float64}}

    ## Reaction network: 
    # X1+X2 → 2X2 (susceptible becomes infected)
    # X2 → X3 (infected becomes recovered)
    ## ME:
    # ∂/∂t{p(X(t)=[x1,x2,x3])}  = κ1⋅(x1+1)⋅(x2-1)⋅p(X(t)=[x1+1,x2-1,x3]) 
    #                           + κ2⋅(x2+1)⋅p(X(t)=[x1,x2+1,x3-1])
    #                           - (κ1⋅x1+κ2)⋅x2⋅p(X(t)=[x1,x2,x3])

    # Gillespie
    x2 = [1]
    x3 = [0]
    x1 = [ntot]
    tvec = [0.0]
    while tvec[end]<tmax
        
        # Sum propensity functions
        a = κ1*x1[end]*x2[end]+κ2*x2[end]

        # Time increment
        Δt = rand(Exponential(1.0/a))

        # Store time
        push!(tvec,tvec[end]+Δt)

        # Which reaction happens
        m = rand()<=(κ1*x1[end]*x2[end]/a) ? 1 : 2

        # Decrease alive counter
        if m==1
            push!(x1,x1[end]-1)
            push!(x2,x2[end]+1)
            push!(x3,x3[end])
        else
            push!(x1,x1[end])
            push!(x2,x2[end]-1)
            push!(x3,x3[end]+1)
        end

    end

    # Return (t,x(t)) tuple
    return (tvec,x1,x2,x3)

end