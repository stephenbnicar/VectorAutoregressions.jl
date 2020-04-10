function coeftable(v::VarEstimate)
    B = coef(v)
    seB = stderror(v)
    m, K = size(B)
    lags = v.lags

    colnms = [
        "Estimate",
        "Std. Error",
        "t value",
    ]

    rownms = Vector{String}()
    for l = 1:lags
        for k = 1:K
            push!(rownms, "y$k.l$l")
        end
    end
    rownms = v.trend ? ["trend"; rownms] : rownms
    rownms = v.constant ? ["Intercept"; rownms] : rownms

    ctable = Array{CoefTable}(undef, K)
    for k = 1:K
        Bk   = B[:, k]
        seBk = seB[:, k]
        tk = Bk ./ seBk
        mat = hcat(Bk, seBk, tk)
        ctable[k] = CoefTable(mat, colnms, rownms)
    end
    return ctable
end
