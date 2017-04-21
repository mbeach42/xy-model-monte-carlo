include("MetropolisXYModelCartesian.jl")

#Compile once
T = Array{Float32}(linspace(0.2, 2, 5))
MetropolisXYModelCartesian.iterate_over_temperatures(T, Int32(5), Int32(6), Int32(100))

function timeit(T::Int32, L::Int32, N::Int32)
  time = zeros(10)
  for i = 1:10
    tic()
    MetropolisXYModelCartesian.iterate_over_temperatures((Array{Float32}(linspace(0.2, 2, T))), (L), Int32(10), (N))
    time[i] = 2*L^2*(10+N)*T/1e6/toq();
  end
  println("Mean time: ", mean(time))
  println("Std. Dev. time: ", sqrt(var(time)))
  println("Max time: ", maximum(time))
  println("Min time: ", minimum(time))
end

timeit(Int32(20), Int32(10), Int32(500))
