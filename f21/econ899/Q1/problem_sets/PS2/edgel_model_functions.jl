

#keyword-enabled structure to hold model primitives
@with_kw struct Primitives
    β::Float64 = 0.9932 #discount rate
    α::Float64 = 1.5 # coefficient of relative risk aversion
    a_min::Float64 = -2 #asset lower bound
    a_max::Float64 = 5 #asset upper bound
    na::Int64 = 1000 #number of asset grid points
    nS::Int64 = 2 #number of states
    a_grid::Array{Float64, 1} = collect(range(a_min, length = na, stop = a_max)) #asset grid
    S::Array{Float64, 1} = collect([1, 0.5]) # earnings in each state
    Π::Array{Float64, 2} = collect([0.97 0.03; 0.05 0.05])
end


#structure that holds model results
mutable struct Results
    val_func::Array{Float64, 2} # value function
    pol_func::Array{Float64, 2} # policy function
    μ::Array{Float64, 2}        # stationary distribution
    q̄::Float64                  # market-clearing bond price
end

#function for initializing model primitives and results
function Initialize()
    prim        = Primitives()                                      # initialize primtiives
    val_func    = zeros(prim.nS, prim.na)                           # initial value function guess
    pol_func    = zeros(prim.nS, prim.na)                           # initial policy function guess
    μ           = repeat([1/(prim.nS*prim.na)], prim.nS, prim.na)   # initial distribution guess
    q̄           = (prim.β + 1)/2                                    # initial price guess
    res         = Results(val_func, pol_func, μ, q̄)                 # initialize results struct
    prim, res #return deliverables
end

# Bellman Operator
function Bellman(prim::Primitives, res::Results)
    @unpack val_func, q̄ = res                   # unpack value function and bond price
    @unpack a_grid, S, β, α, na, nS, Π = prim   # unpack model primitives
    v_next = zeros(nS, na)                      # next guess of value function to fill

    #choice_lower = 1 #for exploiting monotonicity of policy function
    for S_index = 1:nS
        for a_index = 1:na
            a = a_grid[a_index] #value of k
            e = S[S_index] #value of z
            candidate_max = -Inf #bad candidate max
            budget = e + q̄*a #budget

            for ap_index in 1:na #loop over possible selections of k', exploiting monotonicity of policy function
                c = budget - a_grid[ap_index] #consumption given k' selection
                if c>0 #check for positivity
                    exp_val = 0 # expected value in next period
                    for Snext = 1:nS
                        exp_val = exp_val +
                                    Π[S_index, Snext]*val_func[Snext, ap_index]
                    end
                    val = ((c^(1-α)-1)/(1-α)) + β*exp_val #compute value
                    if val>candidate_max #check for new max value
                        candidate_max = val #update max value
                        res.pol_func[S_index, a_index] = a_grid[ap_index] #update policy function
                    end
                end
            end
            v_next[S_index, a_index] = candidate_max #update value function
        end
    end
    v_next #return next guess of value function
end

# Value function iteration
function V_iterate(prim::Primitives, res::Results; tol::Float64 = 1e-4, err::Float64 = 100.0)
    n = 0 #counter

    while err>tol #begin iteration
        v_next = Bellman(prim, res) #spit out new vectors
        err = maximum(abs.(v_next.-res.val_func)) #reset error level
        res.val_func = v_next #update value function
        n+=1
    end
    #println("Value function converged in ", n, " iterations.")
end

# Solve for μ
function Dist_Solve(prim::Primitives, res::Results; tol::Float64 = 1e-4, err::Float64 = 100.0)
    @unpack nS, na, a_grid, Π = prim    # number of states and asset grid points
    @unpack μ, pol_func = res           # policy function and marginal distribution
    n = 0                               # counter

    while err > tol # begin iteration 
        μ₁ = zeros(nS, na)  # next period distribution
        for s = 1:nS # loop through employment states
            for a = 1:na # loop through asset levels for current holdings
                indChoose = findall(pol_func[s, :] .== a_grid[a])
                for sc = 1:nS
                    μ₁[s, a] = μ₁[s, a] + sum(μ[sc, indChoose])*Π[sc, s]
                end
            end
        end
        n+=1
        err = norm(μ₁.-μ) # reset error level

        if (n % 100 == 0)
            #println(n, " iterations; err = ", err)
        end
        μ = μ₁              # Update stationary distribution
    end

    res.μ = μ
    #println("μ converged in ", n, " iterations.")
end

# solve the model
function Solve_model(prim::Primitives, res::Results; tol::Float64 = 1e-3, Ed::Float64 = 100.0, θ::Float64 = 0.8)
    @unpack q̄ = res     # initial price guess
    n = 0               # counter

    while abs(Ed) > tol

        V_iterate(prim, res)    # retrieve policy function for current q̄
        Dist_Solve(prim, res)   # retrieve stationary distribution for current q̄

        Ed = sum(res.μ.*res.pol_func)   # calculate excess demand at current q̄
        if Ed < 0 & (abs(Ed) >= tol)      # adjust price toward bounds according to tuning parameter
            res.q̄ = θ*q̄ + (1-θ)*prim.β
        elseif abs(Ed) >= tol 
            res.q̄ = θ*q̄ + 1 - θ
        end
        n+=1

        if (n % 10 == 0)
            println(n, " iterations; excess demand= = ", Ed)
        end
    end
    println("q̄ converged in ", n, " iterations.")
end
##############################################################################
