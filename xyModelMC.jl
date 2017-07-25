############################################
# XY model Monte Carlo code
# By Matthew Beach
############################################

# This module provides a reasonably efficient Monte Carlo sampling of the XY model in two-dimensions. 
# Guidelines:
# 1) Iterating over temperatures is best done with parallel so be sure to run the .jl script with -p N where N is the number of cores. 
# 2) To avoid saving to a file specify "none" for the directory name.
# 3) 'clusters' is the number of clusters formed in the Wolff algorithm or the number of spin flips/L^2 in the Metropolic-Hasting algorithm

# References:
# The Metropolic-Hasting algorithm used here is specific to the q-state model as discussed here https://journals.aps.org/prb/pdf/10.1103/PhysRevB.26.6201

############################################
# Example:
############################################
# L = 4
# N_eq = 10
# N = 100
# T_low = 0.1
# T_high = 3.0
# num_Ts = 30
# clusters = 10
# dir_name = "none"
# E, M, Cv, X, Stiff = qModelMC.iterate_over_temperatures(T_low, T_high, num_Ts, L, q, N_eq, N, clusters, "wolff", dir_name);
# using PyPlot      # plot the results
# xlabel("Temperature")
# ylabel("Energy")
# plot(linspace(T_low, T_high, num_Ts), E, "." )

@everywhere module xyModelMC

using NPZ

immutable SpinArray
    L::Int                            # lattice size
    config::Array{Float64, 2}           # lattice of q-spin clock states (integers in[0, q-1])
end

magnetization(S::SpinArray) =  abs(sum(sinpi.(2S.config))^2 + sum(cospi.(2S.config))^2 )/S.L^2.

function neighbors(i::Int, j::Int, L::Int) # nearest neighbours
    up = [mod1(i+1,L), j]
    down = [mod1(i-1,L), j]
    left = [i, mod1(j-1,L)]
    right = [i, mod1(j+1,L)]
    return up, down, left, right
end

function energy(S::SpinArray)
    E, M = 0.0, S.config
    for j = 1:S.L, i = 1:S.L
        E -= cospi(2(M[i,j]-M[mod1(i+1,S.L), j])) + cospi(2(M[i,j]-M[i, mod1(j+1,S.L)])) 
    end
    return E/2.0
end

function spin_stiffness_y_direction(S::SpinArray, T::Float64)
    term1, term2 = 0.0, 0.0
    for j = 1:S.L, i = 1:S.L
        x = S.config[i,j] - S.config[i, mod1(j+1,S.L)]
        term1 += cospi(2x)
        term2 += sinpi(2x)
    end
    stiffness =  (term1 - 1./T*term2^2)/(S.L^2)
end

function thermo_quantities(T, L, N_eq, N, clusters, algorithm)
    S = SpinArray(L, rand(L, L)) 
    nrg, mag, stiff, configs = zeros(N), zeros(N), zeros(N), Array{Float64}(N, L^2)
    
    for i = 1:N_eq  # Warm-up
        algorithm == "metro" ? mcsweep!(S, T, clusters) : wolff!(S, T, clusters)
    end
    
    for j = 1:N 
        algorithm == "metro" ? mcsweep!(S, T, clusters) : wolff!(S, T, clusters)
        configs[j,:] = reshape(S.config, L^2)
        nrg[j] = energy(S)
        mag[j] = magnetization(S)
        stiff[j] = spin_stiffness_y_direction(S, T)
    end
    
    E = sum(nrg)/(N*L^2)          
    M = sum(mag)/(N*L^2)                       
    Cv = (dot(nrg,nrg)/N - sum(nrg)^2/N^2)/T^2/L^2
    X = (dot(mag, mag)/N - sum(mag)^2/N^2)/T/L^2
    Stiff = sum(stiff)/N
        
    return E, M, Cv, X, Stiff, configs
end

function iterate_over_temperatures(T_low, T_high, num_Ts, L, N_eq, N, clusters, algorithm, dir_name)
    T_range = linspace(T_low, T_high, num_Ts)
    E, M, Cv, X, Stiff = SharedArray{Float64}(num_Ts),SharedArray{Float64}(num_Ts),SharedArray{Float64}(num_Ts),SharedArray{Float64}(num_Ts),SharedArray{Float64}(num_Ts)
    Configs, configs_array = SharedArray{Float64}(num_Ts*N, L^2), Array{Float64}(num_Ts*N,L^2)

    @sync @parallel for j = 1:num_Ts
         E[j], M[j], Cv[j], X[j], Stiff[j], Configs[(j-1)*N + 1:j*N, :] = thermo_quantities(T_range[j], L, N_eq, N, clusters, algorithm)
    end
    
    if dir_name != "none"
        configs_array = copy(Configs)  # save configs to file
        mkdir(dir_name) 
        mkdir(dir_name*"/L=$L/") 
        npzwrite(dir_name*"/L=$L/"".configs.npy", configs_array)
        save_parameters(T_low, T_high, num_Ts, L, N_eq, N, clusters, algorithm, dir_name*"/L=$L/")
    end

    return E, M, Cv, X, Stiff
end

function save_parameters(T_low, T_high, num_Ts, L, N_eq, N, clusters, algorithm, dir_name)
    open(dir_name*"parameters.txt", "w") do f
        write(f,"L = $L\n")
        write(f,"N_eq = $N_eq\n")
        write(f,"N = $N\n")
        write(f,"clusters = $clusters\n")
        write(f,"algorithm = $algorithm\n")
    end
end

function energy_single_flip(S::SpinArray, trial_spin, i, j)
    E, Etrial, M = 0.0, 0.0, S.config
    for (n,m) in neighbors(i, j, S.L)
        Etrial -= cospi(2(M[i,j]-M[n,m]))
        E      -= cospi(2(trial_spin-M[n,m]))
    end
    return E - Etrial
end

function mcstep!(S::SpinArray, T, x, y)
    old_spin = S.config[x, y]
    trial_spin = rand() # randomly flip one spin left/right one unit
    dE = energy_single_flip(S, trial_spin, x, y)                                              # compute energy of a single spin flip
    if dE < 0 || rand() < exp(-1.0/T*dE)
        S.config[x, y] = trial_spin       # accept the new state
    else
        S.config[x, y] = old_spin         # reject the new state
    end
end

function mcsweep!(S::SpinArray, T, num_metro)
    x, y = rand(1:S.L, 2)
    for i = vcat(x:S.L, 1:x),j = vcat(x:S.L, 1:x), k = 1:num_metro
        mcstep!(S, T, i, j)
    end
end

function wolff!(S::SpinArray, T, num_clusters)
    for iter = 1:num_clusters
        cluster = falses(S.L,S.L)    # create empty cluster
        vectors = zeros(S.L,S.L,2)   # create empty array for spin vectors
        
        v_angle = rand()
        v_vec = [cospi(v_angle), sinpi(v_angle)]
        
        i, j = rand(1:S.L, 2)  # randomly select a lattice point (i,j)
        vectors[i,j,:]  = [cospi(2S.config[i,j]), sinpi(2S.config[i,j])] #compute the sin/cos of that point
        createcluster!(S, T, cluster, i, j, v_vec, vectors)
        flipcluster!(S, cluster, v_vec, vectors)
    end
end

function createcluster!(S::SpinArray, T, cluster, i, j, v_vec, vectors)
    cluster[i,j] = true                                 # assign current spin to cluster
    for (n,m) in neighbors(i, j, S.L)
        if cluster[n,m] == false
            vectors[n,m,:] = [cospi(2S.config[n,m]), sinpi(2S.config[n,m])]           # store spin-vector in array 
            if rand()  < 1. - exp(-2/T*dot(v_vec, vectors[i,j,:])*dot(v_vec, vectors[n,m,:])) 
                createcluster!(S, T, cluster, n, m, v_vec, vectors) 
            end 
        end
    end
end

function flipcluster!(S::SpinArray, cluster, v_vec, vectors)
    shaped_vectors = reshape(vectors, S.L^2, 2)
    for i in find(cluster)
        n_vec = shaped_vectors[i,:]                                   # look up the n-vector
        n_new = n_vec - 2*dot(n_vec, v_vec)*v_vec                     # reflect about v
        S.config[i]  = mod2pi(atan2(n_new[2],n_new[1]))/(2pi)         # compute new spin angle
    end
end

end #end xyModelMC

# compile the model once to test
xyModelMC.iterate_over_temperatures(0.1, 3.0, 2, 2, 1, 2, 1, "wolff", "none");
xyModelMC.iterate_over_temperatures(0.1, 3.0, 2, 2, 1, 2, 1, "metro", "none");

    