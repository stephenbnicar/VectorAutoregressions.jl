using VectorAutoregressions
using CSV, DataFrames
using Test

cd(@__DIR__)

@testset "VAR" begin
    Y = DataFrame!(CSV.File("testdata1.csv"))
    v = VAR(Y, 2)

    @test coef(v) == v.B
    @test stderror(v) == v.seB 
    @test residuals(v) == v.U
    @test fitted(v) == v.Yhat

    Bhat = [-0.017  0.016  0.013
            -0.320  0.044 -0.002
             0.146 -0.153  0.225
             0.961  0.289 -0.264
            -0.161  0.050  0.034
             0.115  0.019  0.355
             0.934 -0.010 -0.022]
    @test isapprox(coef(v), Bhat, atol=2e-3)

    sigu = [0.00212963  7.16167e-5   0.00012324
            7.16167e-5  0.000137338  6.14587e-5
            0.00012324  6.14587e-5   8.92035e-5]
    @test isapprox(v.ΣU, sigu, atol=1e-8)
end
