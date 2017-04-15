Pkg.add("CRlibm")

using CRlibm
CRlibm.setup()  # for CRlibm 4.0

Pkg.clone("https://github.com/musm/Sleef.jl")
using Sleef

function bench(f, r)

    @eval @time Base.$f.($r)

    @eval @time CRlibm.$f.($r, [RoundNearest])

    @eval @time Sleef.$f.($r)
end

function bench(f, r)

    output = similar(r)

    @eval @time $output .= Base.$f.($r)

    @eval @time $output .= CRlibm.$f.($r, [RoundNearest])

    @eval @time $output .= Sleef.$f.($r)
end

function bench(N)
    r = [reinterpret(Float64, rand(Int)) for i in 1:N]

    for f in (:exp, :log, :sin, :cos, :tan, :sinh, :cosh)
        println("$f:")
        if f == :log
            bench(f, r[r .> 0.0])
        else
            bench(f, r)
        end
        println()
    end
end
