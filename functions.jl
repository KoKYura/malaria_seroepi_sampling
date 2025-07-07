using Distributions, DataFrames, GLM, JLD2, Random, Optim, CSV

# Scenario 1
# Generate data sets based on age structures
function get_age_s1(n, p_age1, p_age2, p_age3, p_age4)
    p_age5 = 1.0 - (p_age1 + p_age2 + p_age3 + p_age4)

    n1 = round(Int, n * p_age1)
    n2 = round(Int, n * p_age2)
    n3 = round(Int, n * p_age3)
    n4 = round(Int, n * p_age4)
    n5 = n - (n1 + n2 + n3 + n4)  # ensure total = n

    l_age1 = rand(1:10, n1)
    l_age2 = rand(11:20, n2)
    l_age3 = rand(21:30, n3)
    l_age4 = rand(31:40, n4)
    l_age5 = rand(41:50, n5)

    all_ages = vcat(l_age1, l_age2, l_age3, l_age4, l_age5)

    df_age = DataFrame(ID = 1:n, age = all_ages)
    return df_age
end

# Get age-specific seroprevalence for RCM
function get_sero_positivity_1l(age, l1, rho)
    p = l1 * (1 - exp(-(l1 + rho) * age)) / (l1 + rho)
    return p
end

# Get log-likelihood of RCM estimating lambda only
function get_ll_t1_par1(l1, rho, age, y)
    p = get_sero_positivity_1l.(age, l1, rho)
    return -sum(y .* log.(p) .+ (1 .- y) .* log.(1 .- p))
end

# Get log-likelihood of RCM estimating lambda and rho 
function get_ll_t1_par2(params, age, y)
    l1, rho = params

    penalty = 0.0
    if l1 < 0 || l1 > 1 || rho < 0 || rho > 0.1
        penalty += 1e10
    end
    
    ϵ = 1e-10  
    p = get_sero_positivity_1l.(age, l1, rho)
    p = clamp.(p, ϵ, 1 - ϵ)
    ll = -sum(y .* log.(p) .+ (1 .- y) .* log.(1 .- p))

    return ll + penalty
end

function get_res_s1_par1(n, p_age1, p_age2, p_age3, p_age4, lambda, rho)

    # Data simulation
    df = get_age_s1(n, p_age1, p_age2, p_age3, p_age4)
    df.y_sp = [rand(Binomial(1, get_sero_positivity_1l(age, lambda, rho))) for age in df.age]

    # Parameter estimation by MLE
    ## estimate only lambda
    obj(x) = get_ll_t1_par1(x, rho, df.age, df.y_sp)
    res_par1   = optimize(obj, 0.0, 10.0, Brent())
    est_l = Optim.minimizer(res_par1)

    df_res = DataFrame(
        est_l = est_l
    )
    return df_res
end

function get_res_s1_par2(n, p_age1, p_age2, p_age3, p_age4, lambda, rho)

    # Data simulation
    df = get_age_s1(n, p_age1, p_age2, p_age3, p_age4)
    df.y_sp = [rand(Binomial(1, get_sero_positivity_1l(age, lambda, rho))) for age in df.age]

    # Parameter estimation by MLE
    ## estimate lambda and rho
    res_par2 = optimize(p -> get_ll_t1_par2(p, df.age, df.y_sp),
                           [0.1, 0.1], NelderMead())
    est_l, est_r = Optim.minimizer(res_par2)

    df_res = DataFrame(
        est_l = est_l
    )

    return df_res
end

# Scenario 2
function get_age_s2(n, tau, p_age1, p_age2, p_age3, p_age4)
    p_age5 = 1.0 - (p_age1 + p_age2 + p_age3 + p_age4)

    n1 = round(Int, n * p_age1)
    n2 = round(Int, n * p_age2)
    n3 = round(Int, n * p_age3)
    n4 = round(Int, n * p_age4)
    n5 = n - (n1 + n2 + n3 + n4)  

    l_age1 = rand(1:tau, n1)
    l_age2 = rand(tau+1:tau+10, n2)
    l_age3 = rand(tau+11:tau+20, n3)
    l_age4 = rand(tau+21:tau+30, n4)
    l_age5 = rand(tau+31:tau+40, n5)

    all_ages = vcat(l_age1, l_age2, l_age3, l_age4, l_age5)

    df_age = DataFrame(ID = 1:n, age = all_ages)
    return df_age
end

# Get age-specific sero prevalence for RCM-scSCR
function get_sero_positivity_2l(age::Real, l1, l2, rho, tau)
    if age <= tau
        return get_sero_positivity_1l(age, l2, rho)
    else
        return get_sero_positivity_1l(tau, l2, rho) +
               get_sero_positivity_1l(age - tau, l1, rho) * exp(-(l2 + rho) * tau)
    end
end

# Get log likelihood of RCM-scSCR estimating lambdas only

function get_ll_t2_par2(params, rho, age, y, tau)
    l1, l2 = params

    penalty = 0.0
    if l1 < 0 || l1 > 1 || l2 < 0 || l2 > 1
        penalty += 1e10
    end
    
    ϵ = 1e-10  
    p = get_sero_positivity_2l.(age, l1, l2, rho, tau)
    p = clamp.(p, ϵ, 1 - ϵ)
    ll = -sum(y .* log.(p) + (1 .- y) .* log.(1 .- p))

    return ll + penalty
end

# Get log likelihood of RCM-scSCR estimating lambdas only

function get_ll_t2_par3(params, age, y, tau)
    l1, l2, rho = params

    penalty = 0.0
    if l1 < 0 || l1 > 1 || l2 < 0 || l2 > 1 || rho < 0 || rho > 0.1
        penalty += 1e10
    end
    
    ϵ = 1e-10  
    p = get_sero_positivity_2l.(age, l1, l2, rho, tau)
    p = clamp.(p, ϵ, 1 - ϵ)
    ll = -sum(y .* log.(p) .+ (1 .- y) .* log.(1 .- p))

    return ll + penalty
end

function get_res_s2_par2(n, p_age1, p_age2, p_age3, p_age4, l1, l2, rho, tau)

    # Simulate dataset
    df = get_age_s2(n, tau, p_age1, p_age2, p_age3, p_age4)
    df.y_sp = [rand(Binomial(1, get_sero_positivity_2l(age, l1, l2, rho, tau))) for age in df.age]
    
    # MLE using RCM-scSCR with known SRR
    res_t2_par2 = optimize(p -> get_ll_t2_par2(p, rho, df.age, df.y_sp, tau),
                           [0.1, 0.1], NelderMead())
    est_l1_2param, est_l2_2param = Optim.minimizer(res_t2_par2)

    df_res = DataFrame(
        tau = tau,
        est_l1 = est_l1_2param,
        est_l2 = est_l2_2param,
        )

    return df_res
end

function get_res_s2_par3(n, p_age1, p_age2, p_age3, p_age4, l1, l2, rho, tau)

    # Simulate dataset
    df = get_age_s2(n, tau, p_age1, p_age2, p_age3, p_age4)
    df.y_sp = [rand(Binomial(1, get_sero_positivity_2l(age, l1, l2, rho, tau))) for age in df.age]

    # MLE using RCM-scSCR with unknown SRR
    res_t2_par3 = optimize(p -> get_ll_t2_par3(p, df.age, df.y_sp, tau),
                           [0.1, 0.1, 0.1], NelderMead())
    est_l1_3param, est_l2_3param, est_rho_3param = Optim.minimizer(res_t2_par3)

     df_res = DataFrame(
        tau = tau,
        est_l1 = est_l1_3param,
        est_l2 = est_l2_3param,
        )

    return df_res
end
































