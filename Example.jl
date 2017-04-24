include("MetropolisXYModel.jl")

include("MetropolisXYModelCartesian2.jl")
include("XYPlots.jl")

 L = 8
 N_eq = 200
 N_mc = 2000

T = Array(linspace(0.5, 1.9, 50))
E, Cv, Stiff = MetropolisXYModelCartesian2.iterate_over_temperatures(T, L, N_eq, N_mc)
XYPlots.allplots(T, L, E, Cv, Stiff)
