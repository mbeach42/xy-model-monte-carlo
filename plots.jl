using PyPlot

@everywhere function all_plots_XY(T_range, E, Cv, X)
    fig = figure("my_plot",figsize=(10,4))

    subplot(221)
    xlabel("Temperature")
    ylabel("Energy")
    plot(T_range, E, ".")

    subplot(223)
    xlabel("Temperature")
    ylabel("Cv")
    plot(T_range, Cv, "." ,color = "orange")

    subplot(224)
    xlabel("Temperature")
    ylabel("Sus")
    plot(T_range, X, ".", color = "green");
end

@everywhere function E_M_plot(T_range, E, M)
    fig8 = figure("my_plot",figsize=(10,4))

    subplot(221)
    xlabel("Temperature")
    ylabel("Energy")
    plot(T_range, E, ".", color = "blue")

    subplot(222)
    xlabel("Temperature")
    ylabel("Magnetization")
    plot(T_range, abs(M), ".", color = "red");
end

@everywhere function conv_plot(E, E2)
    fig4 = figure("energy_plot",figsize=(10,8))

    subplot(221)
    xlabel("Monte Carlo time")
    ylabel("Energy of current state")
    plot(E, ".", color = "blue")

    subplot(222)
    xlabel("Monte Carlo time")
    ylabel("Average energy")
    plot(E2, ".", color = "purple");

end

@everywhere function plot_X(T_range, X)
    fig4 = figure("energy_plot",figsize=(10,4))
    xlabel("Temperature")
    ylabel(L"$\chi$")
    plot(T_range, X, ".", color = "orange");
end

@everywhere function plot_Cv(T_range, X)
    fig12 = figure("energy_plot",figsize=(10,4))
    xlabel("Temperature")
    ylabel(L"$C_v$")
    plot(T_range, X, ".", color = "green");
end


###XY plots
@everywhere function plots_XY(T_range,L, E, Cv, Stiff)
    fig6 = figure("my_plot",figsize=(10,5))

    subplot(221)
    xlabel("Temperature")
    ylabel("Energy")
    plot(T_range, E, ".", label = L)
    legend()

    subplot(222)
    xlabel("Temperature")
    ylabel(L"C_v")
    plot(T_range, Cv, ".", label = L)
    #legend()

    subplot(223)
    xlabel("Temperature")
    ylabel(L"$x$-Stiffness")
    #xlim(0.4, 2.0)
    ylim(0,1)
    plot(T_range, T_range*2/pi, color = "black", "-", label = L"\frac{2}{\pi}T")
    plot(T_range, Stiff, ".", label = L);
    #legend()
end

@everywhere function arrow_plot(conf)
    figure("arrow_plot", figsize = (5, 3))
    X = collect(1:size(conf)[1])
    Y = collect(1:size(conf)[1])
    U = sin(conf)
    V = cos(conf)
    xlim(0., size(conf)[1] +1)
    ylim(0., size(conf)[1] +1)
    quiver(X,Y,U,V, cmap = ColorMap("hot")) ;
end
