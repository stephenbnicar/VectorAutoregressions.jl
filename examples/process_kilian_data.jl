using CSV
using Statistics
using DataFrames

cd(@__DIR__)
fedfunds = CSV.read("./source_data/fedfunds.txt"; header=0)
deflator = CSV.read("./source_data/gnpdeflator.txt", header=0)
rgdp = CSV.read("./source_data/realgnp.txt", header=0)

irate = Vector{Float64}()
for i=1:3:length(fedfunds[:,3])
  push!(irate, mean(fedfunds[i:i+2,3]))
end

infl = diff(log.(deflator[:,3]))*100
drgdp = diff(log.(rgdp[:,3]))*100

kilian_data = DataFrame(drgdp = drgdp, irate = irate, infl = infl)

CSV.write("kilian_data.csv", kilian_data)
