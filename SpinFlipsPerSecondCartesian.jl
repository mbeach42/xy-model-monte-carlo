<<<<<<< Updated upstream:SpinFlipsPerSecondCartesian.jl
include("MetropolisXYModelCartesian.jl")
=======
include("MC_2D_XY.jl")

#Compile one for speed
iterate_over_temperatures_XY(Array(linspace(0.6, 0.8, 5)), 5, 5, 5)
#iterate_over_temperatures_XY(linspace(0.6, 0.8, 5), 5, 5, 15)

>>>>>>> Stashed changes:timing-functions.jl

#Compile once
T = Array{Float64}(linspace(0.2, 2, 5))
MetropolisXYModelCartesian.iterate_over_temperatures(T, (5), (6), (100))

function timeit(T, L, N)
  time = zeros(10)
  for i = 1:10
    tic()
<<<<<<< Updated upstream:SpinFlipsPerSecondCartesian.jl
    MetropolisXYModelCartesian.iterate_over_temperatures(Array(linspace(0.2, 2, T)), (L), (10), (N))
=======
    iterate_over_temperatures_XY(Array(linspace(0.2, 2, T)), L, 10, N)
>>>>>>> Stashed changes:timing-functions.jl
    time[i] = 2*L^2*(10+N)*T/1e6/toq();
  end
  println("Mean time: ", mean(time))
  println("Std. Dev. time: ", sqrt(var(time)))
  println("Max time: ", maximum(time))
  println("Min time: ", minimum(time))
end

<<<<<<< Updated upstream:SpinFlipsPerSecondCartesian.jl
timeit((20), (10), (500))
=======
timeit(20, 15, 1000)

#timeit(20, 10, 1000)
>>>>>>> Stashed changes:timing-functions.jl
