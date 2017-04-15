include("MC_2D_XY.jl")
#include("old_XY.jl")
#Compile
iterate_over_temperatures_XY(linspace(0.6, 0.8, 5), 5, 5, 5)


function timeit(T, L, N)
  time = zeros(10)
  for i = 1:10
    #Time
    tic()
    iterate_over_temperatures_XY(linspace(0.2, 2, T), L, 10, N)
    time[i] = 2*L^2*(10+N)*T/1e6/toq();

  end

  println("Mean time: ", mean(time))
  println("Std. Dev. time: ", sqrt(var(time)))
  println("Max time: ", maximum(time))
  println("Min time: ", minimum(time))
end
