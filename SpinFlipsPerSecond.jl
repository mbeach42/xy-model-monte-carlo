include("MetropolisXYModel.jl")

#Compile once
T = Array(linspace(0.2, 2, 5))
MetropolisXYModel.iterate_over_temperatures(T, (5), (6), (100))

function timeit(T::Int, L::Int, N::Int)
  time = zeros(10)
  for i = 1:10
    tic()
    MetropolisXYModel.iterate_over_temperatures((Array(linspace(0.2, 2, T))), (L), (10), (N))
    time[i] = 2*L^2*(10+N)*T/1e6/toq();
  end
  println("Mean time: ", mean(time))
  println("Std. Dev. time: ", sqrt(var(time)))
  println("Max time: ", maximum(time))
  println("Min time: ", minimum(time))
end

timeit((20), (10), (500))
