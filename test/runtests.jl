using BilevelOptimization
using BilevelOptimization.BilevelFlowProblems

using Test
import LinearAlgebra

using Cbc: CbcSolver
using JuMP

test_bp() = BilevelLP(
    [1.,0.], [0.],
    [-.1 0.], -ones(1,1), [-1.],
    [-1.],
    [-.1 0.], ones(1,1), ones(1,)
)

@testset "BilevelLP type" begin
    bp = test_bp()
    @test bp isa BilevelLP{V,M} where{V<:AbstractVector{Float64},M<:AbstractMatrix{Float64}}
end

@testset "test basic problem" begin
    bp = test_bp()
    (m, x, y, λ, _) = build_blp_model(bp, CbcSolver())
    status = JuMP.solve(m)
    @test status === :Optimal
    xv = JuMP.getvalue(x)
    yv = JuMP.getvalue(y)
    @test xv[1] ≈ 0.0
    @test yv[1] ≈ 1.0
end

@testset "Integrality is registered" begin
    bp = test_bp()
    push!(bp.Jx, 2) # second useless variable is integer
    (m, x, y, λ, _) = build_blp_model(bp, CbcSolver())
    @test JuMP.getcategory(x[2]) == :Int
    @test JuMP.getcategory(x[1]) == :Cont
    @test JuMP.getcategory(y[1]) == :Cont
end

@testset "Problem from Dempe-Mersha 2006" begin
    # problem taken from:
    # A.G. Mersha, S. Dempe,
    # Linear bilevel programming with upper-level constraints depending on lower-level solution.
    # Appl. Math. Computation (180), 2006
    cx = [-1.]
    cy = [-2.]
    G = zeros(2,1) .+ [-2.,1.]
    H = zeros(2,1) .+ [3.,1.]
    q = [12.,14.]
    d = [-1.]
    A = zeros(2,1) .+ [-3.,3.]
    B = ones(2,1)
    b = [-3.,30.]
    prob = BilevelLP(
        cx, cy,
        G, H, q,
        d, A, B, b
    )
    (m, x, y, λ, _) = build_blp_model(prob, CbcSolver())
    st = JuMP.solve(m)
    @test st === :Optimal
    @test getvalue(x)[1] ≈ 8.
    @test getvalue(y)[1] ≈ 6.

    # problem modification: passing upper-level constraints to lower level
    A = zeros(4,1) .+ [-3.,3.,-2.,1.]
    B = ones(4,1); B[3,1] = 3.
    b = [-3.,30.,12.,14.]
    prob = BilevelLP(
        cx, cy,
        G, H, q,
        d, A, B, b
    )
    (m, x, y, λ, _) = build_blp_model(prob, CbcSolver())
    st = solve(m)
    @test st === :Optimal
    @test getvalue(x)[1] ≈ 6.
    @test getvalue(y)[1] ≈ 8.
end

@testset "Problem with upper bounds" begin
    # problem taken from:
    # Bilevel Programming Problems
    # S. Dempe, V. Kalashnikov, G. Pérez-Valdés, N. Kalashnykova
    # Springer 2015
    cx = [2.,1.]
    cy = [2.,-1.]
    G = zeros(0,2)
    H = zeros(0,2)
    q = Float64[]
    d = [0.,0.]
    F = [1. 0.;0. 1.]
    A = zeros(1,2)
    B = [-2. 1.]
    b = [0.]
    prob = BilevelLP(
        cx, cy,
        G, H, q,
        d, A, B, b, Int[], F
    )
    setlowerbound(prob, BilevelOptimization.upper, 1, -1.)
    setupperbound(prob, BilevelOptimization.upper, 1, 1.)
    setlowerbound(prob, BilevelOptimization.upper, 2, -1.)
    setupperbound(prob, BilevelOptimization.upper, 2, -0.75)
    setlowerbound(prob, BilevelOptimization.lower, 1, -Inf64)
    setupperbound(prob, BilevelOptimization.lower, 1, 2.)
    setlowerbound(prob, BilevelOptimization.lower, 2, 0.)
    setupperbound(prob, BilevelOptimization.lower, 2, 2.)
    @test size(prob.B) == (prob.ml,prob.nl)
    (m, x, y, λ, _) = build_blp_model(prob, CbcSolver())
    st = JuMP.solve(m)
    @test st === :Optimal
    @test all(getvalue(x) .≈ (-1.,-1.))
end

@testset "Upper-level integer problem" begin
    # problem taken from:
    # Bilevel Programming Problems
    # S. Dempe, V. Kalashnikov, G. Pérez-Valdés, N. Kalashnykova
    # Springer 2015
    cx = [2.,0.]
    cy = [3.,2.,6.]
    G = [ 4.  1.
         -4. -1.]
    H = zeros(2,3)
    q = [10.,-10.]
    Jx = Int64[1, 2]
    d = [-5.,-8.,-1.]
    F = zeros(2,3)
    B = [4. 2. 0.
         2. 4. 1.
        -1. 0. 0.
         0. -1. 0.
         0. 0. -1.
         ]
    A = [-1.  0.
          0. -1.
          0.  0.
          0.  0.
          0.  0.]
    b = [0.,0.,0.,0.,0.]
    intprob = BilevelLP(
        cx, cy,
        G, H, q,
        d, A, B, b, Jx, ylowerbound = false
    )
    (m, x, y, λ, _) = build_blp_model(intprob, CbcSolver())
    st = JuMP.solve(m)
    @test st === :Optimal
    @test all(getvalue(x) .≈ (2.,2.))
    @test all(getvalue(y) .≈ (1//3,1//3,0))
    B = [4. 2. 0.
         2. 4. 1.]
    A = [-1.  0.
          0. -1.]
    b = [0.,0.]
    intprob = BilevelLP(
        cx, cy,
        G, H, q,
        d, A, B, b, Jx
    )
    (m, x, y, λ, _) = build_blp_model(intprob, CbcSolver())
    st = JuMP.solve(m)
    @test st === :Optimal
    @test all(getvalue(x) .≈ (2.,2.))
    @test all(getvalue(y) .≈ (1//3,1//3,0))
end

function test_bflow()
    init_cost = [
        0. 1. 1. 4.
        0. 0. 0. 1.
        0. 0. 0. 1.
        0. 0. 0. 0.
    ]
    taxable_edges = [
        false true true false
        false false false true
        false false false true
        false false false false
    ]
    tax_options = zeros(4,4,5)
    for i in 1:4, j in 1:4
        if taxable_edges[i,j]
            tax_options[i,j,:] .= (0.0, 0.5, 1.0, 1.5, 2.0)
        end
    end

    capacities = [
        0. 3. 2. 3.
        0. 0. 0. 2.
        0. 0. 0. 1.
        0. 0. 0. 0.
    ]
    BilevelFlowProblem(init_cost,taxable_edges,capacities,tax_options,3.)
end

@testset "Flow problem matrix construction" begin
    bfp = test_bflow()
    (B,b) = BilevelFlowProblems.flow_constraint_standard(bfp)
    @test size(B)[1] == bfp.nv-2 + 1 + bfp.nv*bfp.nv
    @test sum(B[1,:]) ≈ -3.
    @test B[1,5] ≈ -1.
    @test B[1,9] ≈ -1.
    @test B[1,13] ≈ -1.
    @test b[1] ≈ -bfp.minflow
    # flow conservation
    @test B[2,5] ≈ 1.
    @test B[2,14] ≈ -1.
    @test b[2] ≈ 0.
    @test B[3,9] ≈ 1.
    @test B[3,15] ≈ -1.
    @test b[3] ≈ 0.
    # capacity constraints
    @test all(LinearAlgebra.I - B[4:end,:] .≈ 0)
    @test all(b[4:end] .≈ [bfp.capacities[i] for i in eachindex(bfp.capacities)])
end

@testset "Bilevel flow JuMP model" begin
    bfp = test_bflow()
    (m, r, y, f, λ) = build_blp_model(bfp, CbcSolver())
    st = JuMP.solve(m)
    @test st === :Optimal
    @test getobjectivevalue(m) ≈ 6.
    for j in 1:size(r)[2]
        for i in 1:size(r)[1]
            @test getvalue(r[i,j]) ≈ sum(getvalue(y[i,j,:]).*bfp.tax_options[i,j,:]) * getvalue(f[i,j])
        end
    end
end

@testset "Complementarity types" begin
    m = JuMP.Model()
    # testing non-crashing
    BilevelOptimization.add_complementarity_constraint(m, SOS1Complementarity, [], [], [], [])
    BilevelOptimization.add_complementarity_constraint(m, SOS1Complementarity, [], [], [], [])
end

@testset "Basic problem with specified methods" begin
    bp = test_bp()
    (m, x, y, λ) = build_blp_model(bp, CbcSolver(), comp_method = SOS1Complementarity)
    status = JuMP.solve(m)
    @test status === :Optimal
    xv = JuMP.getvalue(x)
    yv = JuMP.getvalue(y)
    @test xv[1] ≈ 0.0
    @test yv[1] ≈ 1.0
    (m, _, _, _) = build_blp_model(bp, CbcSolver(), comp_method = SOS1Complementarity())
    status = JuMP.solve(m)
    @test status === :Optimal
    # arbitrary big-enough bounds
    (m, x, y, λ) = build_blp_model(bp, CbcSolver(), comp_method = BoundComplementarity(100.,100.))
    status = JuMP.solve(m)
    @test status === :Optimal
    xv = JuMP.getvalue(x)
    yv = JuMP.getvalue(y)
    @test xv[1] ≈ 0.0
    @test yv[1] ≈ 1.0
end

@testset "Bilevel flow big-M bounds" begin
    bfp = test_bflow()
    (m, r, y, f, λ) = build_blp_model(bfp, CbcSolver(), comp_method = BoundComplementarity(100., 100.))
    st = JuMP.solve(m)
    @test st === :Optimal
    @test getobjectivevalue(m) ≈ 6.
    for j in 1:size(r)[2]
        for i in 1:size(r)[1]
            @test getvalue(r[i,j]) ≈ sum(getvalue(y[i,j,:]).*bfp.tax_options[i,j,:]) * getvalue(f[i,j])
        end
    end
end

@testset "Bilevel flow big-M vector bounds" begin
    bfp = test_bflow()
    bounds_method = BoundComplementarity(10 .* ones(19), 30.)
    (m, r, y, f, λ) = build_blp_model(bfp, CbcSolver(), comp_method = bounds_method)
    st = JuMP.solve(m)
    @test st === :Optimal
    @test getobjectivevalue(m) ≈ 6.
    for j in 1:size(r)[2]
        for i in 1:size(r)[1]
            @test getvalue(r[i,j]) ≈ sum(getvalue(y[i,j,:]).*bfp.tax_options[i,j,:]) * getvalue(f[i,j])
        end
    end
    bounds_method = BoundComplementarity(10., 10 .* ones(19))
    (m, r, y, f, λ) = build_blp_model(bfp, CbcSolver(), comp_method = bounds_method)
    st = JuMP.solve(m)
    @test st === :Optimal
    @test getobjectivevalue(m) ≈ 6.
    for j in 1:size(r)[2]
        for i in 1:size(r)[1]
            @test getvalue(r[i,j]) ≈ sum(getvalue(y[i,j,:]).*bfp.tax_options[i,j,:]) * getvalue(f[i,j])
        end
    end
end

@testset "Bilevel flow big-M infeasible" begin
    bfp = test_bflow()
    (m, r, y, f, λ) = build_blp_model(bfp, CbcSolver(), comp_method = BoundComplementarity(0.1, 0.1))
    st = JuMP.solve(m)
    @test st === :Infeasible
end
