include("MC_2D_XY.jl")
include("plots.jl")
include("timing-functions.jl")

T = Array(linspace(0.1, 2.5, 50))
E,M,S = iterate_over_temperatures_XY(T, 3, 5, 5)
E,M,S = iterate_over_temperatures_XY(T, 6, 500, 1000)
all_plots_XY(T, E,M,S)



println("\n # of T = 20, L = 10, N = 500 ")
timeit(20, 10, 500)

println("\n # of T = 20, L = 10, N = 1000 ")
timeit(20, 10, 1000)
