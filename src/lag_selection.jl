struct LagSelectionCriteria
    table::DataFrame
    selection::Dict
end

function lagselect()

end

function show(io::IO, ls::LagSelectionCriteria)
    println(io, "VAR Lag Selection Criteria")
    println(io, "---------------------------")
    println(io, "Criteria   Lag")
    println(io, "---------------------------")
    for k in keys(ls.selection)
        println(io, rpad(string(k), 11), ls.selection[k])
    end
    println(io, "---------------------------")
end
