function coeftable(v::VarEstimate)
    B = coef(v)
    seB = stderror(v)
    m, K = size(B)
    lags = v.lags
    varnames = v.varnames
    obs = v.obs
    dofr = obs - m

    colnms = [
        "Estimate",
        "Std. Error",
        "t value",
        "Pr > |t|",
    ]

    rownms = Vector{String}()
    for l = 1:lags
        for k = 1:K
            push!(rownms, "$(varnames[k]).l$l")
        end
    end
    rownms = v.trend ? ["trend"; rownms] : rownms
    rownms = v.constant ? ["Intercept"; rownms] : rownms

    ctable = Array{CoefTable}(undef, K)
    for k = 1:K
        Bk   = round.(B[:, k], sigdigits=4)
        seBk = round.(seB[:, k], sigdigits=4)
        tk = Bk ./ seBk
        pk = 2 * ccdf.(TDist(dofr), abs.(tk))
        # mat = hcat(Bk, seBk, tk)
        mat = hcat(Bk, seBk, tk, pk)
        # ctable[k] = CoefTable(mat, colnms, rownms)
        ctable[k] = CoefTable(mat, colnms, rownms, 4, 3)
    end
    return ctable
end
