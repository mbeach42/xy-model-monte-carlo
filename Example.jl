include("MetropolisXYModel.jl")
include("XYPlots.jl")

const L = 6
const N_eq = 100
const N_mc = 5000

T = Array(linspace(0.1, 2.5, 50))
E, Cv, Stiff = MetropolisXYModel.iterate_over_temperatures(T, L, N_eq, N_mc)
XYPlots.allplots(T, L, E, Cv, Stiff)
