β = coef(var_df)[:, 1]
σ = stderror(var_df)[:, 1]
t = β ./ σ
# p = 2 * ccdf.(TDist(dof_residual(obj)), abs.(t))
mat = hcat(β, σ, t)
# lims = (100 * (1 - level) / 2, 100 * (1 - (1 - level) / 2))
colnms = [
    "PE   ",
    "SE   ",
    "t-value",
]
rownms = ["x$i" for i in 1:length(β)]
CoefTable(mat, colnms, rownms)
