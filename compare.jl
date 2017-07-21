include("MetropolisXYModel.jl")
include("MetropolisXYModelCartesian.jl")

L = 6
T = 0.89
N_eq = 0
N_steps = 50


config = MetropolisXYModel.random_config(L)
X, Y = cospi(config), sinpi(config)
E1, stiff1  = zeros(N_steps), zeros(N_steps)
E2, stiff2  = zeros(N_steps), zeros(N_steps)
x, y = rand(1:L, 2*L^2, N_steps+N_eq), rand(1:L, 2*L^2, N_steps+N_eq)
rs, deltaTheta = rand(2*L^2, N_steps+N_eq), 0.25*pi*rand(2*L^2, N_steps+N_eq)

E1[1] = MetropolisXYModel.energy(config, L)
E2[1] = MetropolisXYModelCartesian.energy(X, Y, L)
a,b = 0,0
for i = 1:N_steps
    config = MetropolisXYModel.random_config(L)
    config2 = deepcopy(config)
    X,Y = cospi(config), sinpi(config)
    x,y = rand(1:L, N_steps), rand(1:L, N_steps)
    rs = 2*rand(N_steps)
    deltaTheta = rand(N_steps)*pi

    for i = 1:20
        f = rand(1:N_steps)
        MetropolisXYModel.mcstep!(config, T, L, x[f],y[f],rs[f],deltaTheta[f])
        MetropolisXYModelCartesian.mcstep!(config2, X, Y, T, L, x[f],y[f],rs[f],deltaTheta[f])
        if maximum(abs(config-config2)) > 1e-10
            println("maximun: ", maximum(abs(config-config2)) )
            println("Spherical E: ", MetropolisXYModel.energy_single_flip(config, deltaTheta[f],x[f],y[f],L))
            E = MetropolisXYModelCartesian.energy_single_flip(X, Y, x[f], y[f], L) #compute energy of spin flip
            dx = cospi(config2[x[f],y[f]] + deltaTheta[f]) - X[x[f],y[f]]
            dy = sinpi(config2[x[f],y[f]] + deltaTheta[f]) - Y[x[f],y[f]]
            config2[x[f],y[f]] += deltaTheta[f]
            X[x[f],y[f]] += dx #randomly flip one spin
            Y[x[f],y[f]] += dy
            Eflip = MetropolisXYModelCartesian.energy_single_flip(X, Y, x[f], y[f], L) #compute energy of spin flip

            deltaE = Eflip - E
            println("Cartesian E: ", deltaE)
            println(x[f],y[f],rs[f],deltaTheta[f])
        end
    end

        #MetropolisXYModel.mcstep!(config2, T, L, x,y,rs,deltaTheta)
        #MetropolisXYModelCartesian.mcstep!(config, X, Y, T, L, x,y,rs,deltaTheta)

    #println(maximum(cospi(config2) - X))
    E1[i] = MetropolisXYModel.energy(config, L)
    stiff1[i] = MetropolisXYModel.stiffness_x(config, T, L)

    E2[i] = MetropolisXYModelCartesian.energy(X, Y, L)
    stiff2[i] = MetropolisXYModelCartesian.stiffness_x(X, Y, T, L)
end

#display(E1-E2)

println("energy diff: ", maximum(abs(E1-E2)))
