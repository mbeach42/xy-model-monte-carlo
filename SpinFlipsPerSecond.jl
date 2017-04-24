include("MetropolisXYModel.jl")

#Compile one for speed
MetropolisXYModel.iterate_over_temperatures(Array(linspace(0.6, 0.8, 5)), 5, 5, 5)

function timeit(T, L, N)
  time = zeros(10)
  Ts = Array(linspace(0.2, 2.0, T))
  println(Ts)

  for i = 1:10
    tic()
    MetropolisXYModel.iterate_over_temperatures(Ts, L, 10, N)
    time[i] = 2*L^2*(10+N)*T/1e6/toq();
  end

  println("Mean time: ", mean(time))
  println("Std. Dev. time: ", sqrt(var(time)))
  println("Max time: ", maximum(time))
  println("Min time: ", minimum(time))
end


timeit(20, 10, 500)
