module XYPlots

using PyPlot

export allplots, arrow_plot

function allplots(T_range, L, E, Cv, Stiff)
    fig = figure("my_plot", figsize=(10,5))

    subplot(221)
    xlabel("Temperature")
    ylabel("Energy")
    plot(T_range, E, ".", color = "red", label = label = "L = $(L)")
    legend()

    subplot(222)
    xlabel("Temperature")
    ylabel("C_v")
    plot(T_range, Cv, ".")

    subplot(223)
    xlabel("Temperature")
    ylabel("x-Stiffness")
    ylim(0,1)
    plot(T_range, T_range*2/pi, color = "black", "-", label = "\frac{2}{\pi}T")
    plot(T_range, Stiff, ".");

    return fig
end

function arrowplot(conf)
    fig = figure("arrow_plot", figsize = (5, 3))
    X = collect(1:size(conf)[1])
    Y = collect(1:size(conf)[1])
    U = sin(conf)
    V = cos(conf)
    xlim(0., size(conf)[1] +1)
    ylim(0., size(conf)[1] +1)
    quiver(X,Y,U,V, cmap = ColorMap("hot")) ;
    return fig
end

end #XYPlots
