@everywhere function periodic_index(a::Int64, L::Int64)::Int64 #This impliments periodic boudnary conditions (PBC's)
    if a == 0
        return L
    elseif a == L + 1
        return 1
    else
      return a
  end
end

@fastmath @everywhere function random_config_XY(L::Int64)::Array{Float64,2} #generate a random (L,L) matrix of +/-1
    @inbounds return 2.0*pi*rand(L,L)
end

@fastmath @everywhere function energy_singleflip_XY(M::Array{Float64,2}, deltaTheta::Float64, i::Int64, j::Int64, L::Int64)::Float64 #energy of configuration
    E, Etrial = 0.0, 0.0
    @fastmath @inbounds Etrial += - cos(M[i,j]-M[i, periodic_index(j+1,L)]) - cos(M[i,j]-M[i, periodic_index(j-1, L)]) - cos(M[i,j]-M[periodic_index(i+1, L),j]) - cos(M[i,j]-M[periodic_index(i-1, L),j])
    @inbounds M[i,j] -= deltaTheta
    @fastmath @inbounds E = - cos(M[i,j]-M[i, periodic_index(j+1,L)]) - cos(M[i,j]-M[i, periodic_index(j-1, L)]) - cos(M[i,j]-M[periodic_index(i+1, L),j]) - cos(M[i,j]-M[periodic_index(i-1, L),j])
    @inbounds M[i,j] += deltaTheta
    return Etrial - E
end

@fastmath @everywhere function energy_XY(M::Array{Float64,2}, L::Int64)::Float64 #energy of configuration
    E = 0.0
    @inbounds for j = 1:L, i = 1:L #Sum the energy of the four nearest-neighbour spins
        @inbounds E += - cos(M[i,j]-M[i, periodic_index(j+1,L)])
        @inbounds E += - cos(M[i,j]-M[i, periodic_index(j-1, L)])
        @inbounds E += - cos(M[i,j]-M[periodic_index(i+1, L),j])
        @inbounds E += - cos(M[i,j]-M[periodic_index(i-1, L),j])
    end
    return E/4.0
end

@fastmath @everywhere function stiffness_XY(M::Array{Float64,2}, T::Float64, L::Int64)::Float64 # in y-direction
    term1, term2 = 0.0, 0.0
    @inbounds for j = 1:L, i = 1:L
        @inbounds term1 += cos(M[i,j]-M[i, periodic_index(j+1,L)])
        @inbounds term2 += sin(M[i,j]-M[i, periodic_index(j+1,L)])
    end
    stiffness =  term1 - 1./T*term2^2
    return stiffness *= 1./(L^2)
end

@fastmath @everywhere function MC_step_XY!(config::Array{Float64,2}, T::Float64, L::Int64,
                                x::Int64, y::Int64, r::Float64, deltaTheta::Float64)#Monte Carlo step
    config[x, y] += deltaTheta #randomly flip one spin
    deltaE = energy_singleflip_XY(config, deltaTheta, x, y, L) #compute energy of spin flip
    if deltaE > 0
      return nothing
    else
      config[x, y] = mod(config[x, y], 2.0*pi)
      r > exp(-1.0/T*deltaE) ? config[x, y] -= deltaTheta : nothing
    end
end

@fastmath @everywhere function MC_sweep_XY!(config::Array{Float64,2}, T::Float64, L::Int64, x::Array{Int64, 1}, y::Array{Int64, 1}, rs::Array{Float64, 1}, deltaTheta::Array{Float64, 1}) #Sweep over L^2 MC steps
    for i = 1:2*L^2
        MC_step_XY!(config, T, L, x[i], y[i], rs[i], deltaTheta[i])
    end
end

@fastmath @everywhere function thermo_quantities_XY(T::Float64, L::Int64, N_eq::Int64, N_steps::Int64)::Tuple{Float64,Float64,Float64}

    config = random_config_XY(L)
    E, stiff  = zeros(N_steps), zeros(N_steps)

    x, y= rand(1:L, 2*L^2, N_steps+N_eq), rand(1:L, 2*L^2, N_steps+N_eq)
    rs, deltaTheta = rand(2*L^2, N_steps+N_eq), 0.25*pi*rand(2*L^2, N_steps+N_eq) + 2.0*atan(T)

    for i = 1:N_eq #Run a few times to equilibriate
         MC_sweep_XY!(config, T, L, x[:,i], y[:,i], rs[:,i], deltaTheta[:,i])
    end

    for i = 1:N_steps #Runs which will be included in the thermodynamic averages
        MC_sweep_XY!(config, T, L, x[:,N_eq+i], y[:,N_eq+i], rs[:,N_eq+i], deltaTheta[:,N_eq+i])
        E[i] = energy_XY(config, L)
        stiff[i] = stiffness_XY(config, T, L)
    end

    E_avg = sum(E)/(N_steps*L^2) #average energy at T
    Cv = (dot(E,E)/N_steps - sum(E)^2/N_steps^2)/(T^2)
    Stiff = sum(stiff)/N_steps
    return E_avg, Cv, Stiff
end


@fastmath @everywhere function iterate_over_temperatures_XY(Ts::LinSpace{Float64}, L::Int64,
    N_eq::Int64, N_steps::Int64)::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}}

    F = length(Ts)
    E, Cv, Stiff = zeros(F), zeros(F), zeros(F) #SharedArray(Float64,F), SharedArray(Float64,F), SharedArray(Float64,F)

    for j = 1:length(Ts)
        E[j], Cv[j], Stiff[j] = thermo_quantities_XY(Ts[j], L, N_eq, N_steps)
    end
    return E, Cv, Stiff
end
