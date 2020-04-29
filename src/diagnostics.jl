"""
    sic(v::VarEstimate)

Return the Schwarz (Bayesian) Information Criterion for VAR model `v`.
"""
function sic(v::VarEstimate)
    obs, K = size(residuals(v))
    nparam = size(coef(v), 1)
    ΣU = v.ΣU
    ΣUml = ((obs - nparam) / obs) * ΣU
    return logdet(ΣUml) + (log(obs) / obs) * nparam * K
end

"""
    hqc(v::VarEstimate)

Return the Hannan-Quinn Criterion for VAR model `v`.
"""
function hqc(v::VarEstimate)
    obs, K = size(residuals(v))
    nparam = size(coef(v), 1)
    ΣU = v.ΣU
    ΣUml = ((obs - nparam) / obs) * ΣU
    return logdet(ΣUml) + (2 * log(log(obs)) / obs) * nparam * K
end

struct StabilityCheck
    varnames::Array{String}
    lags::Int
    isstable::Bool
    eigenvals::Array{Number}
    eigenmod::Array{Float64}
end

function checkstable(v::VarEstimate)
    # See Lutkepohl (2006) p.15ff
    p = v.lags
    B = coef(v)
    nparam, K = size(B)
    startidx = nparam - K*p + 1
    A = B[startidx:end, :]
    Acomp = comp_matrix(A, p)
    eigenvals = eigvals(Acomp)
    eigenmod = abs.(eigenvals)
    isstable = all(eigenmod .< 1)
    StabilityCheck(v.varnames, v.lags, isstable, eigenvals, eigenmod)
end

function show(io::IO, stab::StabilityCheck)
    K = length(stab.varnames)
    E = stab.eigenvals
    Emod = stab.eigenmod
    println(io, "VAR Stability Check")
    println(io, "===================")
    print(io, "Endogenous variables: ")
    for k = 1:K-1
        print(io, "$(stab.varnames[k]), ")
    end
    println(io, "$(stab.varnames[K])")
    println(io, "Lags: $(stab.lags)")
    if stab.isstable
        println(io, "VAR is stable")
    else
        println(io, "VAR is not stable")
    end
    println(io, "--------------------------")
    println(io, " Eigenvalues      Modulus")
    println(io, "--------------------------")
    for i in 1:length(E)
        println(io, lpad(string(round(E[i]; digits=3)), 16), lpad(string(round(Emod[i]; digits=3)), 9))
    end
    println(io, "--------------------------")
end
