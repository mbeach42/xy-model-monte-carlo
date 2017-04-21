include("MetropolisXYModel.jl")
include("MetropolisXYModelCartesian.jl")

L = 6
a = 2*rand(L, L)
x, y = cospi(a), sinpi(a)

d = 1.76552

dy = sinpi(a[3,3]+d) - y[3,3]
dx = cospi(a[3,3]+d) - x[3,3]

#?? x = (x^2-y^2-x)/(x^2-y^2-2*x+1), dy = -y/(x^2-y^2-2*x+1)


sin(theta+Dtheta) = sin(theta)cos(DTHETA) + cos(theta)snsin(DTHETA)
= y*cos(dtheta) + x*sin(dtheta) = sin(theta+dtheta)


a[3,3]
a[3,3] += d

x[3,3]
x[3,3] += dx

y[3,3]
y[3,3] += dy

E1, E2 = MetropolisXYModel.energy_single_flip(a, d, 3,3, L)
E3, E4 = MetropolisXYModelCartesian.energy_single_flip(x, y, dx, dy, 3,3, L)

E1 - E3

E2 - E4
############################
a[3,3]

x[3,3]

y[3,3]
#####################
a[3,3] -= d

x[3,3] -= dx

y[3,3] -= dy
