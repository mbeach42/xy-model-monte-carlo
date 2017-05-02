addprocs(6)
println(nprocs(), " nodes")

println("Importing...")

include("MetropolisXYModel.jl")
include("XYPlots.jl")

println("Run once for compiling...")
name2 = "/home/matt/github/xy-model-monte-carlo/xy_data_L=$(4)_N_eq=$(10)_N=$(10)"
isdir(name2) ? nothing : mkdir(name2)
iterate_over_temperatures_parallel(Array(linspace(0.1, 4.0, 10)), 4, 10, 10);

println("Declaring...")
L = 32
N_eq = 2000
N = 10000
println(L, " ", N_eq, " ", N)
name1 = "/home/matt/github/xy-model-monte-carlo/xy_data_L=$(L)_N_eq=$(N_eq)_N=$(N)"
isdir(name1) ? nothing : mkdir(name1)


println("Starting...")

tic()
T = Array(linspace(0.2, 2.2, 41))
E, Cv, Stiff = iterate_over_temperatures_parallel(T, L, N_eq, N);

t = toq()


allplots(T, L, E, Cv, Stiff)


println("Time is ", 2*L^2*(N_eq+N)*length(T)/1e6/t)
println("Time for L = 32, N = 10,000 is ", t*32^2/L^2*100000/N/60/60, " hours.")
