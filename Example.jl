include("MetropolisXYModelCartesian2.jl")
#include("XYPlots.jl")

 L = 6
 N_eq = 100
 N_mc = 5000

T = Array(linspace(0.1, 2.5, 50))
E, Cv, Stiff = MetropolisXYModelCartesian2.iterate_over_temperatures(T, L, N_eq, N_mc)
#XYPlots.allplots(T, L, E, Cv, Stiff)
