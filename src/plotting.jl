## Plot recipe
@recipe function f(irf::IRF; shock = 1, response = 1, ci = true)
    point = irf.point[resp, :, shock]
    upper = irf.upper[resp, :, shock]
    lower = irf.lower[resp, :, shock]
    h = length(point)
    x = 0:1:(h-1)
    if ci
        ymax = maximum(upper)
        ymin = minimum(lower)
    else
        ymax = maximum(point)
        ymin = minimum(point)
    end
    cush = max(abs(0.1*ymin), abs(0.1*ymax))

    ylims  := (ymin-cush, ymax+cush)
    xlims  := (0, h-0.75)
    xticks := 0:1:h-1
    legend := false
    linewidth --> 2
    linecolor --> :black

    ## Try substituting hline! call here instead
    @series begin
        linewidth := 1
        linecolor := :darkgray
        x, zeros(h)
    end

    @series begin
        linestyle := :solid
        x, point
    end

    if ci
        @series begin
            linewidth := 1
            linestyle := :dash
            x, [upper lower]
        end
    end

end
