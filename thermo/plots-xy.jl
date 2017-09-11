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

function spin_stiffness_x_direction(S::SpinArray, T::Float64)
    term1, term2 = 0.0, 0.0
    for j = 1:S.L, i = 1:S.L
        x = S.config[i,j] - S.config[i, mod1(j+1,S.L)]
        term1 += cospi(2x)
        term2 += sinpi(2x)
    end
    stiffness =  (term1 - 1./T*term2^2)/(S.L^2)
end

function loaded_thermo_quantities(T, L, N, config)
    nrg, mag, stiff = zeros(N), zeros(N), zeros(N)
    configs = reshape(config, N,L,L)
    for j = 1:N 
        S = SpinArray(L, configs[j,:,:])
        nrg[j] = energy(S)
        mag[j] = magnetization(S)
        stiff[j] = spin_stiffness_x_direction(S, T)
    end
    E = sum(nrg)/(N*L^2)          
    M = sum(mag)/(N*L^2)                       
    Cv = (dot(nrg,nrg)/N - sum(nrg)^2/N^2)/T^2/L^2
    X = (dot(mag, mag)/N - sum(mag)^2/N^2)/T/L^2
    Stiff = sum(stiff)/N
        
    return E, M, Cv, X, Stiff
end

function loaded_iterate_over_temperatures(T_low, T_high, num_Ts, L, N, configs)
    T_range = linspace(T_low, T_high, num_Ts)
    E, M, Cv, X, Stiff = zeros(num_Ts),zeros(num_Ts),zeros(num_Ts),zeros(num_Ts),zeros(num_Ts)
    for j = 1:num_Ts
        E[j], M[j], Cv[j], X[j], Stiff[j] = loaded_thermo_quantities(T_range[j], L, N, configs[(j-1)*N + 1:j*N, :])
    end
    return E, M, Cv, X, Stiff
end


N = 1000
T_low = 0.1
T_high = 3.0
num_Ts = 64



for L in [8,16,24,32,40,48,56,64,72]
        println("\nL = $L")
        dir_name = "/home/mbeach/xy-paper-1/plotting-xy/"
        ispath(dir_name) == true ? nothing : mkdir(dir_name) 

        configs = npzread("/scratch/mbeach/10-seeds-xy-1000/L=$L/configs.npy")
        E, M, Cv, X, Stiff = loaded_iterate_over_temperatures(T_low, T_high, num_Ts, L, N, configs);

        T_range = linspace(T_low, T_high, num_Ts)
        writedlm(dir_name*"T_$L.txt", T_range, ",")
        writedlm(dir_name*"E_$L.txt", E, ",")
        writedlm(dir_name*"M_$L.txt", M, ",")
        writedlm(dir_name*"Cv_$L.txt", Cv, ",")
        writedlm(dir_name*"X_$L.txt", X, ",")
        writedlm(dir_name*"stiff_$L.txt", Stiff, ",")
end


