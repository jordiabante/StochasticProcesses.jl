#################################################################################
# Dependencies
#################################################################################
using Plots
using StochasticProcesses

#################################################################################
## Pure birth process w/ constant intensity/propensity.
## Equivalent to Poisson process: X(t)~Po(λt).
#################################################################################
# Configuration
λ = 1e1
tmax = 100.0

# Average trajectory: E[X(t)]=λt
p = plot(1:tmax,λ.*collect(1:tmax),lw=4,color="black");

# Plot trajectories of pure birth process
for _ in 1:20
    t,x = StochasticProcesses.sim_traj_birth_ct(λ,tmax)
    plot!(t,x);
end
plot!(xlab="t",ylab="X(t)",title="Constant birth process (λ=$(λ))",legend=false);
display(p)
savefig("/Users/jordiabante/Desktop/Constant-Birth-Process.pdf")

#################################################################################
## Pure birth process w/ linear intensity/propensity.
## Equivalent to Binomial process: X(t)~ TODO
#################################################################################
# Configuration
λ = 5e-1
tmax = 10.0

# Average trajectory: E[X(t)]=exp{λt}
p = plot(1:0.1:tmax,exp.(λ.*collect(1:0.1:tmax)),lw=4,color="black");

# Plot trajectories of pure birth process
for _ in 1:20
    t,x = StochasticProcesses.sim_traj_birth_lin(λ,tmax)
    plot!(t,x);
end
plot!(xlab="t",ylab="X(t)",title="Linear birth process (λ=$(λ))",legend=false);
display(p)
savefig("/Users/jordiabante/Desktop/Linear-Birth-Process.pdf")

#################################################################################
## Pure death process w/ constant intensity/propensity.
## Equivalent to TODO
#################################################################################
# Configuration
μ = 1e-1
x0 = 100
tmax = 1000.0

# Average trajectory: E[X(t)]=max{0,X0-μt}
p = plot(0:tmax,x0.-μ.*collect(0:tmax),lw=4,color="black")

# Plot trajectories of pure birth process
for _ in 1:10
    t,x = StochasticProcesses.sim_traj_death_ct(μ,x0)
    plot!(t,x)
end
plot!(xlab="t",ylab="X(t)",title="Constant death process (μ=$(μ))",legend=false)
display(p)

#################################################################################
## Pure death process w/ linear intensity/propensity.
## Equivalent to Binomial process: X(t)~Binomial(exp{-μt};X0).
#################################################################################
# Configuration
μ = 1e-1
x0 = 100
tmax = 100.0

# Average trajectory: E[X(t)]=λt
p = plot(1:tmax,x0*exp.(-μ.*collect(1:tmax)),lw=4,color="black")

# Plot trajectories of pure birth process
for _ in 1:20
    t,x = StochasticProcesses.sim_traj_death_lin(μ,x0)
    plot!(t,x)
end
plot!(xlab="t",ylab="X(t)",title="Linear death process (μ=$(μ))",legend=false)
display(p)

#################################################################################
## Pure death process w/ linear intensity/propensity.
## Equivalent to Binomial process: X(t)~Binomial(exp{-μt};X0).
#################################################################################
# Configuration
x0 = 1
λ = 1e0
μ = 5e-1
tmax = 10.0

# Average trajectory: E[X(t)]=exp{(λ-μ)t}
p = plot(0:0.1:tmax,x0*exp.((λ-μ).*collect(0:0.1:tmax)),lw=4,color="black");

# Plot trajectories of pure birth process
for _ in 1:20
    t,x = StochasticProcesses.sim_traj_birth_death_lin(λ,μ,tmax)
    plot!(t,x);
end
plot!(xlab="t",ylab="X(t)",title="Linear birth-death process (λ=$(λ),μ=$(μ))",legend=false);
display(p)
savefig("/Users/jordiabante/Desktop/Linear-Birth-Death-Process.pdf")
