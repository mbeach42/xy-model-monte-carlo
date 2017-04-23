include("MetropolisXYModel.jl")
include("MetropolisXYModelCartesian2.jl")

L = 6
a = 2*rand(L, L)
x, y = cospi(a), sinpi(a)

d = 1.76552

dy = sinpi(a[3,3]+d) - y[3,3]
dx = cospi(a[3,3]+d) - x[3,3]

sinpi(d)





a[3,3]
a[3,3] += d

x[3,3]
x[3,3] += dx

y[3,3]
y[3,3] += dy

E0 = MetropolisXYModel.energy_single_flip(a, d, 3,3, L)

Eflip = MetropolisXYModelCartesian.energy_single_flip(x, y, 3,3, L)
x[3,3] -= dx
y[3,3] -= dy
Eflip = MetropolisXYModelCartesian.energy_single_flip(x, y, 3,3, L)

Eflip - E
E0

#
#
####
