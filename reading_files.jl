L = 10
N_eq = 0
N = 5000
temps = Array(linspace(0.2, 2.2, 21))

path = "/home/matt/github/xy-model-monte-carlo/xy_data_L=10_N_eq=500_N=5000/"
files = readdir(path)
#config = readdlm(path*files[1])

Stiff = zeros(length(temps))
stiff = zeros(N)

function periodic_index(a::Int64, L::Int64)::Int64 #This impliments periodic boudnary conditions (PBC's)
   if a == 0
       return L
   elseif a == L + 1
       return 1
   else
     return a
  end
end


function stiffness_x(M,T::Float64,L::Int64)::Float64 # in y-direction
   term1, term2 = 0.0, 0.0
    for j = 1:L, i = 1:L
       x = M[i,j] - M[periodic_index(i+1,L), j]
        term1 += cospi(x)
        term2 += sinpi(x)
   end
   stiffness =  term1 - 1./T*term2^2
   return stiffness *= 1./(L^2)
end


for j = 1:length(temps)
    T = temps[j]
    println("T is ", T)
    println(files[j])
    config = readdlm(path*files[j])
    #display(config)
    for i = 1:N
        M = reshape(config[i,:],L,L)
        #display(M)
        stiff[i] = stiffness_x(M, T, L)
        #config = vcat(readdlm(path*files[i]), config)
    end
    Stiff[j] = sum(stiff)/N
end

using PyPlot
plot(temps, Stiffls
  ls
  )
