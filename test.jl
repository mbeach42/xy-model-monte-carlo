#include("MC_2D_XY.jl")
include("timing-functions.jl")
include("old_XY.jl")
include("plots.jl")

T = linspace(0.5, 2.5,50)
L = 5

E, Cv, S = iterate_over_temperatures_XY(T, L, 500, 2000)

all_plots_XY(T, E, Cv, S)
