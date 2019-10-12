cd(@__DIR__)
using CSV
using Statistics
using DataFrames

fedfunds = CSV.read("fedfunds.txt"; header=0)
deflator = CSV.read("gnpdeflator.txt", header=0)
rgdp = CSV.read("realgnp.txt", header=0)

irate = Vector{Float64}()
for i=1:3:length(fedfunds[:,3])
  push!(irate, mean(fedfunds[i:i+2,3]))
  # irate=[irate; mean(fedfunds[i:i+2,3])]
end

infl = diff(log.(deflator[:,3]))*100
drgdp = diff(log.(rgdp[:,3]))*100

kilian_data = DataFrame(drgdp = drgdp, irate = irate, infl = infl)

CSV.write("kilian_data.csv", kilian_data)
