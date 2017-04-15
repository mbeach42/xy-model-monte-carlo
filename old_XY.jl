using PyPlot
using JLD
@everywhere function periodic_index(a::Int64, L::Int64)::Int64 #This impliments periodic boudnary conditions (PBC's)
    if a == 0
        return L
    elseif a == L + 1
        return 1
    else
      return a
  end
end


@everywhere function random_config_XY(L::Int64)::Array{Float64,2} #generate a random (L,L) matrix of +/-1
    return 2.0*pi*rand(L,L)
end

@everywhere function energy_singleflip_XY(M::Array{Float64,2}, a::Float64, i::Int64, j::Int64, L::Int64)::Float64 #energy of configuration
    E1, E2 = 0.0, 0.0
    #E trial
    E2 = -cos(M[i,j]-M[i, periodic_index(j+1, L)]) -cos(M[i,j]-M[i, periodic_index(j-1, L)])-cos(M[i,j]-M[periodic_index(i+1, L),j]) -cos(M[i,j]-M[periodic_index(i-1, L),j])

    #E original
    M[i,j] -= a
    E1 = -cos(M[i,j]-M[i, periodic_index(j+1, L)]) -cos(M[i,j]-M[i, periodic_index(j-1, L)])-cos(M[i,j]-M[periodic_index(i+1, L),j]) -cos(M[i,j]-M[periodic_index(i-1, L),j])

    M[i,j] += a
    return E2 - E1
end

@everywhere function energy_XY(M::Array{Float64,2}, L::Int64)::Float64 #energy of configuration
    E = 0.0
    for i = 1:L
        for j = 1:L #Sum the energy of the four nearest-neighbour spins
            E += -cos(M[i,j]-M[i, periodic_index(j+1, L)]) -cos(M[i,j]-M[i, periodic_index(j-1, L)])-cos(M[i,j]-M[periodic_index(i+1, L),j]) -cos(M[i,j]-M[periodic_index(i-1, L),j])
        end
    end
    return E/4.0
end

@everywhere function stiffness_XY(M::Array{Float64,2}, T::Float64, L::Int64)::Float64 # in y-direction
    term1, term2 = 0.0, 0.0
    for i = 1:L
        for j = 1:L
            term1 += cos(M[i,j]-M[i, periodic_index(j+1,L)])
            term2 += sin(M[i,j]-M[i, periodic_index(j+1,L)])
        end
    end
    stiffness =  term1 - 1./T*term2^2
    stiffness *= 1./(L^2)
end

@everywhere function MC_step_XY!(config::Array{Float64,2}, T::Float64, L::Int64)::Array{Float64, 2} #Monte Carlo step
    flipx, flipy = rand(1:L), rand(1:L)
    #a = rand()*2.0*pi
    #a = mod(0.25*pi*randn() + 2.0*atan(T), 2.0*pi)
    a = 0.25*pi*randn() + 2.0*atan(T)
    config[flipx, flipy] += a #randomly flip one spin
    deltaE = 1.0*energy_singleflip_XY(config, a, flipx, flipy, L) #compute energy of spin flip
    if rand() < min(1.0, exp(-1.0/T*deltaE)) #if spin flip is favourable, accept new state
        return config
    else
        config[flipx, flipy] -= a #otherwise reverse the spin flip
        return config
    end
end

@everywhere function MC_sweep_XY!(config::Array{Float64,2}, T::Float64, L::Int64)::Array{Float64, 2} #Sweep over L^2 MC steps
    for i = 1:2*L^2
        config = MC_step_XY!(config, T,  L)
    end
    return config
end

@everywhere function thermo_quantities_XY(T::Float64, L::Int64, N_eq::Int64,N_steps::Int64)::Tuple{Float64,Float64,Float64}

    config = random_config_XY(L)
    E, stiff  = zeros(N_steps), zeros(N_steps)

    for i = 1:N_eq #Run a few times to equilibriate
        config = MC_sweep_XY!(config, T, L)
    end

    for i = 1:N_steps #Runs which will be included in the thermodynamic averages
        config = MC_sweep_XY!(config, T, L)
        E[i] = energy_XY(config, L)
        stiff[i] = stiffness_XY(config, T, L)
    end

    E_avg = sum(E)/(N_steps*L^2) #average energy at T
    Cv = (dot(E,E)/N_steps - sum(E)^2/N_steps^2)/(T^2)
    Stiff = sum(stiff)/N_steps
    return E_avg, Cv, Stiff
end


@everywhere function iterate_over_temperatures_XY(Ts::LinSpace{Float64}, L::Int64,
    N_eq::Int64, N_steps::Int64)::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}}

    F = length(Ts)
    E, M, Stiff = SharedArray(Float64,F), SharedArray(Float64,F), SharedArray(Float64,F)

    @sync @parallel for j = 1:length(Ts)
        E[j], M[j], Stiff[j] = thermo_quantities_XY(Ts[j], L, N_eq, N_steps)
    end
    return E, M, Stiff
end
