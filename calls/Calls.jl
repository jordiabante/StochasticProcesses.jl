#################################################################################
# Dependencies
#################################################################################
using Plots
using StatsPlots
using StochasticProcesses

#################################################################################
## Pure birth process w/ constant intensity/propensity.
## Equivalent to Poisson process: X(t)~Po(λt).
#################################################################################
# Configuration
λ = 1e-1
tmax = 100.0

# Average trajectory: E[X(t)]=λt
p = plot(1:tmax,λ.*collect(1:tmax),lw=4,color="black");

# Plot trajectories of pure birth process
num_trans = []
total_time = []
for _ in 1:20
    t,x = StochasticProcesses.sim_traj_birth_ct(λ,tmax)
    push!(num_trans,length(x))
    push!(total_time,t[end])
    plot!(t,x,seriestype=:step,label="");
end
plot!(xlab="t",ylab="X(t)",title="Constant birth process (λ=$(λ))",legend=false);
display(p)
savefig("/Users/jordiabante/Desktop/Constant-Birth-Process.pdf")

# MLE of λ
λhat = sum(num_trans)/sum(total_time)

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
## Birt and death process - closed system
#################################################################################
# Configuration
κ1 = 1e0
κ2 = 1e0
ntot = 1000
tmax = 4.0

# Average trajectory
p1 = plot(0:0.1:tmax,ntot.*(κ2/(κ1+κ2)).*(1.0.-exp.(-(κ1+κ2).*collect(0:0.1:tmax))),
    ylim=(0,ntot),ylab="Population",color="black",label="E[X1]",lwd=5,
    title="Closed birth-death process (κ1=$(κ1),κ2=$(κ2),N=$(ntot))");
p2 = plot(0:0.1:tmax,ntot.-ntot.*(κ2/(κ1+κ2)).*(1.0.-exp.(-(κ1+κ2).*collect(0:0.1:tmax))),
    ylim=(0,ntot),xlab="Time (t)",ylab="Population",color="black",label="E[X2]",lwd=5);

# Plot trajectories of pure birth process
for _ in 1:20
    t,x2 = StochasticProcesses.sim_traj_birth_death_lin_closed(κ1,κ2,ntot,tmax)
    plot!(p1,t,ntot.-x2,alpha=0.25,label="",seriestype=:step);
    plot!(p2,t,x2,alpha=0.25,label="",seriestype=:step);
end
plot(p1,p2,layout=(2,1),size=(700,600))
savefig("/Users/jordiabante/Desktop/Closed-Birth-Death-Process-$(ntot).pdf")

#################################################################################
## Birth and death process - open system
#################################################################################
# Configuration
κ1 = 8e1
κ2 = 1e0
tmax = 5.0

# Average trajectory
p = plot(0:0.1:tmax,(κ1/κ2).*(1.0.-exp.(-κ2.*collect(0:0.1:tmax))),
    ylim=(0,1.3*κ1/κ2),ylab="Population",color="black",label="E[X]",lwd=5,
    xlab="Time (t)",title="Open birth-death process (κ1=$(κ1),κ2=$(κ2))");

# Plot trajectories of pure birth process
for _ in 1:20
    t,x = StochasticProcesses.sim_traj_birth_death_lin_open(κ1,κ2,tmax)
    plot!(p,t,x,alpha=0.25,label="",seriestype=:step);
end
plot(p,size=(700,400))
savefig("/Users/jordiabante/Desktop/Open-Birth-Death-Process.pdf")

#################################################################################
## SIR model
#################################################################################
# Configuration
κ1 = 2.5e-1
κ2 = 7.5e-1
ntot = 100
tmax = 5.0

# Skeleton plot
p1 = plot(xlim=(0,tmax),ylab="Population",xlab="Days",title="SIR (κ1=$(κ1),κ2=$(κ2),N=$(ntot))");
plot!(p1,[0.0],[0.0],alpha=0.25,label="Susceptible",color="blue");
plot!(p1,[0.0],[0.0],alpha=0.25,label="Infected",color="red");
plot!(p1,[0.0],[0.0],alpha=0.25,label="Recovered",color="green");

# Plot trajectories of pure birth process
for _ in 1:20
    t,x1,x2,x3 = StochasticProcesses.sim_traj_sir(κ1,κ2,ntot,tmax)
    plot!(p1,t,x1,alpha=0.25,label="",color="blue",seriestype=:step);
    plot!(p1,t,x2,alpha=0.25,label="",color="red",seriestype=:step);
    plot!(p1,t,x3,alpha=0.25,label="",color="green",seriestype=:step);
end
# plot(p1,size=(700,400))
# savefig("/Users/jordiabante/Desktop/SIR.pdf")

# Get data
Ms = collect(1)
κ1_mat = fill(NaN,(length(Ms),100))
κ2_mat = fill(NaN,(length(Ms),100))
for i=1:length(Ms)
    for j in 1:100
        
        # Generate trajectories
        data = []
        for _ in 1:Ms[i]
            push!(data,StochasticProcesses.sim_traj_sir(κ1,κ2,ntot,tmax))
        end

        # Perform statistical inference
        κ1hat,κ2hat = StochasticProcesses.mle_sir(data)

        # Fill matrix
        κ1_mat[i,j] = κ1hat
        κ2_mat[i,j] = κ2hat

    end
end

# Plot
κ1_vec = vcat(κ1_mat'...)
κ2_vec = vcat(κ2_mat'...)
κ_vec = vcat(κ1_vec,κ2_vec)
M_vec = vcat([fill(j,100) for j in Ms]...)
M_vec = vcat(M_vec,M_vec)
p_vec = vcat(fill("Estimated κ1",length(κ1_vec)),fill("Estimated κ2",length(κ2_vec)))
xlabel = "Number of observed trajectories (M)"
ylabel = "Parameter estimates"
p2 = hline([κ1],color=:blue,lwd=4,label="True κ1",
    title="MLE in SIR (κ1=$(κ1); κ2=$(κ2); N=$(ntot); T=$(tmax))");
hline!([κ2],color=:orange,lwd=4,label="True κ2");
boxplot!(M_vec,κ_vec,group=p_vec,color=[:blue :orange],xlab=xlabel,
    ylab=ylabel,xlim=(Ms[1]-1,Ms[end]+1),ylim=(0,1.5*max(κ1,κ2)));
# savefig("/Users/jordiabante/Desktop/MLE-SIR-N-$(ntot)-T-$(tmax).pdf")
plot(p1,p2,size=(700,800),layout=(2,1))
# savefig("/Users/jordiabante/Desktop/SIR-Trajectories-Inference.pdf")
