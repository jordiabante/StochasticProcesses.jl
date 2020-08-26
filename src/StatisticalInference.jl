function mle_sir(data)
    
    # Delta of x's
    Δ1 = [-1,1,0]
    Δ2 = [0,-1,1]

    # Get numerator of κ1
    num_1 = 0
    num_2 = 0
    den_1 = 0.0
    den_2 = 0.0
    for m=1:length(data)
        # Get data from trajectory
        t,x1,x2,x3 = data[m]
        # Get finite indices
        val_ind = t .< Inf
        t = t[val_ind]
        x1 = x1[val_ind]
        x2 = x2[val_ind]
        x3 = x3[val_ind]
        # Get differentials
        Δx = [ [x1[i]-x1[i-1],x2[i]-x2[i-1],x3[i]-x3[i-1]] for i=2:length(t)]
        Δt = [ t[i]-t[i-1] for i=2:length(t)]
        # Numerator of MLE
        num_1 += sum(Δx .== [Δ1])
        num_2 += sum(Δx .== [Δ2])
        # Denominator
        den_1 += sum(Δt .* x1[1:(end-1)] .* x2[1:(end-1)] )
        den_2 += sum(Δt .* x2[1:(end-1)] )
    end
    
    # Return  MLE
    num_1/den_1,num_2/den_2

end