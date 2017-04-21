module MetropolisXYModelCartesian

export iterate_over_temperatures, energy, energy_single_flip, periodic_index


function periodic_index(a::Int64, L::Int64)::Int64 #This impliments periodic boudnary conditions (PBC's)
   if a == 0
       return L
   elseif a == L + 1
       return 1
   else
     return a
  end
end

 function random_config(L::Int64)::Tuple{Array{Float64,2},Array{Float64,2}} #generate a random (L,L) matrix of +/-1
    theta = 2.0*rand(L,L)
    return cospi(theta), sinpi(theta)
end

 @fastmath function energy_single_flip(X, Y, deltaX,deltaY, i, j, L)#energy of configuration
   E, Etrial = 0.0, 0.0
   Etrial -= X[i,j]*X[i, periodic_index(j+1,L)] + Y[i,j]*Y[i, periodic_index(j+1,L)]
   Etrial -= X[i,j]*X[i, periodic_index(j-1,L)] + Y[i,j]*Y[i, periodic_index(j-1,L)]
   Etrial -= X[i,j]*X[periodic_index(i+1,L), j] + Y[i,j]*Y[periodic_index(i+1,L), j]
   Etrial -= X[i,j]*X[periodic_index(i-1,L), j] + Y[i,j]*Y[periodic_index(i-1,L), j]
   X[i,j] -= deltaX
   Y[i,j] -= deltaY
   E -= X[i,j]*X[i, periodic_index(j+1,L)] + Y[i,j]*Y[i, periodic_index(j+1,L)]
   E -= X[i,j]*X[i, periodic_index(j-1,L)] + Y[i,j]*Y[i, periodic_index(j-1,L)]
   E -= X[i,j]*X[periodic_index(i+1,L), j] + Y[i,j]*Y[periodic_index(i+1,L), j]
   E -= X[i,j]*X[periodic_index(i-1,L), j] + Y[i,j]*Y[periodic_index(i-1,L), j]
   X[i,j] += deltaX
   Y[i,j] += deltaY
   return Etrial, E
end

  @fastmath function energy(X::Array{Float64,2}, Y, L::Int64)::Float64 #energy of configuration
   E = 0.0
    @inbounds for j = 1:L, i = 1:L #Sum the energy of the four nearest-neighbour spins
        E -= (X[i,j]-X[i, periodic_index(j+1,L)])
        E -= (X[i,j]-X[i, periodic_index(j-1,L)])
        E -= (X[i,j]-X[periodic_index(i+1,L), j])
        E -= (X[i,j]-X[periodic_index(i-1,L), j])
   end
   return E/4.0
end

  @fastmath function stiffness_x(X::Array{Float64,2},Y::Array{Float64,2}, T::Float64, L::Int64)::Float64 # in y-direction
   term1, term2 = 0.0, 0.0
    for j = 1:L, i = 1:L
        x = X[i,j] - X[i, periodic_index(j+1,L)]
        y = Y[i,j] - Y[i, periodic_index(j+1,L)]
        term1 += x
        term2 += y
   end
   stiffness =  term1 - 1./T*term2^2
   return stiffness *= 1./(L^2)
 end

 function mcstep!(X::Array{Float64,2},Y::Array{Float64,2}, T::Float64, L::Int64,
                  x::Int64, y::Int64, r::Float64, deltaX::Float64, deltaY::Float64)#Monte Carlo step
   X[x, y] += deltaX #randomly flip one spin
   Y[x, y] += deltaY
   deltaE = energy_single_flip(X, Y, deltaX, deltaY, x, y, L) #compute energy of spin flip
   if deltaE < 0 #accept new lower energy
      return nothing
   else
     r > exp(-1.0/T*deltaE) ? X[x, y] -= deltaX : nothing
     r > exp(-1.0/T*deltaE) ? Y[x, y] -= deltaY : nothing
   end
end

 function mcsweep!(X::Array{Float64,2},Y::Array{Float64,2}, T::Float64, L::Int64, x::Array{Int64, 1}, y::Array{Int64, 1}, rs::Array{Float64, 1}, deltaX::Array{Float64, 1}, deltaY::Array{Float64, 1}) #Sweep over L^2 MC steps
   for i = 1:2*L^2
       mcstep!(X, Y, T, L, x[i], y[i], rs[i], deltaX[i], deltaY[i])
   end
end

 function thermo_quantities(T::Float64, L::Int64, N_eq::Int64, N_steps::Int64)#::Tuple{Float64,Float64,Float64}
   X, Y = random_config(L)
   E, stiff  = zeros(N_steps), zeros(N_steps)
   x, y = rand(1:L, 2*L^2, N_steps+N_eq), rand(1:L, 2*L^2, N_steps+N_eq)
   rs, deltaTheta = rand(2*L^2, N_steps+N_eq), 0.25*rand(2*L^2, N_steps+N_eq)
   deltaX, deltaY = cospi(deltaTheta), sinpi(deltaTheta)
   @inbounds for i = 1:N_eq #Run a few times to equilibriate
        mcsweep!(X,Y, T, L, x[:,i], y[:,i], rs[:,i], deltaX[:,i], deltaY[:,i])
   end
   @inbounds for i = 1:N_steps #Runs which will be included in the thermodynamic averages
       mcsweep!(X,Y, T, L, x[:,N_eq+i], y[:,N_eq+i], rs[:,N_eq+i], deltaX[:,i], deltaY[:,i])
       E[i] = energy(X, L)
       stiff[i] = stiffness_x(X, Y, T, L)
   end
   E_avg = sum(E)/(N_steps*L^2) #average energy at T
   Cv = (dot(E,E)/N_steps - sum(E)^2/N_steps^2)/(T^2)
   Stiff = sum(stiff)/N_steps
   return E_avg, Cv, Stiff
end


 function iterate_over_temperatures(Ts::Array{Float64, 1}, L::Int64,
   N_eq::Int64, N_steps::Int64)::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}}
   F = length(Ts)
   E, Cv, Stiff = zeros(F), zeros(F), zeros(F) #SharedArray(Float64,F), SharedArray(Float64,F), SharedArray(Float64,F)
   @inbounds for j = 1:length(Ts)
       E[j], Cv[j], Stiff[j] = thermo_quantities(Ts[j], L, N_eq, N_steps)
   end
   return E, Cv, Stiff
end

end
