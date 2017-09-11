println(ARGS[1])
L = parse(Int64, ARGS[1])

include("xyModelMC.jl")

println(L)

N_eq = 500 # WHAT IS THE PROPER EQ TIME FOR XY MODEL WITH CLUSTERS?
N = 1000
T_low = 0.1
T_high = 3.0
num_Ts = 64
clusters = 101

dir_name = "/scratch/mbeach/10-seeds-xy-$N/"
ispath(dir_name) == true ? nothing : mkdir(dir_name) 

@time E, M, Cv, X, Stiff = xyModelMC.iterate_over_temperatures(T_low, T_high, num_Ts, L, N_eq, N, clusters, "wolff", dir_name);


N = 100
dir_name = "/scratch/mbeach/10-seeds-xy-test-$N/"
ispath(dir_name) == true ? nothing : mkdir(dir_name) 

@time E, M, Cv, X, Stiff = xyModelMC.iterate_over_temperatures(T_low, T_high, num_Ts, L, N_eq, N, clusters, "wolff", dir_name);

