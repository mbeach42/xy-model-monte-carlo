include("MetropolisXYModel.jl")
include("MetropolisXYModelCartesian2.jl")

L = 5
T = 0.2
a = 2*rand(L, L)
x, y = cospi(a), sinpi(a)

d = rand()


#E0 = MetropolisXYModel.energy_single_flip(a, d, 3,3, L)
#Eflip = MetropolisXYModel.energy_single_flip(a, d, 3,3, L)

#MetropolisXYModel.energy(a,L)
#MetropolisXYModelCartesian2.energy(x,y,L)
b = copy(a)


sum(cospi(a) - x)
sum(sinpi(a) - y)

sum(a - b)

MetropolisXYModel.mcstep!(a, T, L, 3, 4, 0.765, d)
sum(cospi(a) - x)
sum(cospi(a) - y)
sum(a - b)



MetropolisXYModelCartesian2.mcstep!(b, x, y, T, L, 3, 4, 0.765, d)


println(sum(abs(cospi(a)-x) ))
println(sum(abs(sinpi(a)-y) ))


println(MetropolisXYModel.energy(a, L) - MetropolisXYModelCartesian2.energy(x,y,L))


MetropolisXYModel.thermo_quantities(1.1,5, 100, 1000)

MetropolisXYModelCartesian2.thermo_quantities(1.1,5, 100, 1000)
