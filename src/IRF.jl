struct IRF
    point::Array
    upper::Array
    lower::Array
end

function IRF(v::VarEstimate, h; boot::Bool = false, reps::Int = 500, ci::Float64 = 0.95,
             cumulate = false, nm = false)
    Y = v.data
    p = v.lags
    constant = v.constant
    trend    = v.trend
    B = v.B
    U = v.U
    Σ = v.ΣU

    K = size(B, 2)
    ϕ = simple_irf(B, K, p, h+1)
    Θ = orthogonalize_irf(ϕ, Σ, K, h)
    point = permutedims(Θ, [1, 3, 2])

    if boot
        upper, lower = bootstrap_irf_ci(Y, p, constant, trend, B, U, h; reps=reps, ci=ci)
    else
        upper = zeros(K, h+1, K)
        lower = zeros(K, h+1, K)
    end
    # Normalize responses
    if nm
        for k = 1:K
            upper[:, :, k] = upper[:, :, k]/point[k, 1, k]
            lower[:, :, k] = lower[:, :, k]/point[k, 1, k]
            point[:, :, k] = point[:, :, k]/point[k, 1, k]
        end
    end
    # Cumulative impulse responses
    if cumulate
        upper = cumsum(upper, dims=2)
        lower = cumsum(lower, dims=2)
        point = cumsum(point, dims=2)
    end

    IRF(point, upper, lower)
end

function simple_irf(B, K, p, h)
    nparam = size(B, 1)
    A = B[(nparam-K*p)+1:end, :]

    #   phi(:,:,i) is all simple responses to all shocks at horizon i
    phi = zeros(K, K, h+1)
    # For the calculation, see Lutkepohl, p.52
    bigA = comp_matrix(A, p)
    bigJ = [Matrix{Float64}(I, K, K) zeros(K, K*p-K)]
    for i = 0:h
        phi[:, :, i+1] = bigJ*bigA^i*bigJ'
    end
    return phi
end

function orthogonalize_irf(phi, sigma, K, h)
    S = Symmetric(sigma)
    # Orthogonalize using lower Cholesky factorization
    P = cholesky(S).L
    Θ = zeros(K, K, h+1)
    Θ[:, :, 1] = P

    for i = 1:(h+1)
        Θ[:, :, i] = phi[:, :, i]*P
    end
    return Θ
end

function bootstrap_irf_ci(Y, p, c, t, B, U, h; reps=500, ci=0.95)
    nobs, K  = size(Y)
    uobs = nobs - p
    ps = Y[1:p, :]
    lb = (1-ci)/2
    ub = 1 - lb
    sim_sirf = []
    sim_oirf = []
    sim_irf_store = zeros(K*reps, h+1, K)
    upper = zeros(K, h+1, K)
    lower = zeros(K, h+1, K)

    for rep = 1:reps
        Ysim     = simulate_var(B, K, p, c, t, uobs; presample = ps, resid = U)
        simVAR   = VAR(Ysim, p; constant = c, trend = t)
        simB     = simVAR.B
        simSigma = simVAR.ΣU
        sim_sirf = simple_irf(simB, K, p, h+1)
        sim_oirf = orthogonalize_irf(sim_sirf, simSigma, K, h)
        #=============
        Store results:
        -- "Layer" k (3rd dim) contains all responses to a shock to varible k
        -- Within each layer, all replications for the response of each variable are grouped
        together; e.g. rows [1:reps] contain all replications for the response of variable 1,
        rows [reps+1:2*reps] contain all replications for the response of variable 2, etc.
        ==============#
        for sh = 1:K
            for res = 1:K
                sim_irf_store[reps*(res-1) + rep, :, sh] = sim_oirf[res, sh, :]
            end
        end # storage loop

        # Get quantiles
        for sh = 1:K
            for res = 1:K
                for hor = 1:h+1
                    lower[res, hor, sh] = quantile(sim_irf_store[reps*(res-1)+1:(reps*res), hor, sh], lb)
                    upper[res, hor, sh] = quantile(sim_irf_store[reps*(res-1)+1:(reps*res), hor, sh], ub)
                end
            end
        end # quantile loop
    end # rep loop

    return upper, lower
end

function simulate_var(B, K, p, constant, trend, uobs; presample = nothing,
                      sigma = nothing, resid = nothing)
    # Initialize simulated data
    simY = zeros(uobs+p, K)
    # Set presample values
    if !isa(presample, Nothing)
        simY[1:p, :] = presample
    end
    # Set covariance matrix for disturbances
    if isa(sigma, Nothing)
        Σ = Matrix{Float64}(I, K, K)
    else
        Σ = sigma
    end
    # Draw sample of residuals (with replacement) if provided
    if !isa(resid, Nothing)
        boot = true
        randidx = sample(collect(1:uobs), uobs)
        bootU = resid[randidx, :]
    else
        boot = false
    end

    for t = (p+1):(uobs+p)
        # lagged endogenous variables for RHS
        Zt = transpose(simY[t-1:-1:t-p, :])[:]
        if trend
            Zt = [t; Zt]
        end
        if constant
            Zt = [1; Zt]
        end

        if boot
            simY[t, :] = B'*Zt + bootU[t-p, :]
        else
            simY[t, :] = B'*Zt + rand(MvNormal(zeros(K), Σ))
        end
    end # rep loop

    return simY
end
