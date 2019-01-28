using BilevelOptimization
using BilevelOptimization.BilevelFlowProblems

using Test

using Cbc: CbcSolver
import JuMP

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
    (m, x, y, λ) = build_blp_model(bp, CbcSolver())
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
    (m, x, y, λ) = build_blp_model(bp, CbcSolver())
    @test JuMP.getcategory(x[2]) == :Int
    @test JuMP.getcategory(x[1]) == :Cont
    @test JuMP.getcategory(y[1]) == :Cont
end

function test_bflow()
    init_cost = ones(4, 4)
    init_cost[1,4] = 4.
    taxable_edges = [
        false false true false
        false false false false
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
        0. 0. 0. 0. # these should not be used
    ]
    BilevelFlowProblem(init_cost,taxable_edges,capacities,tax_options,3.)
end

@testset "Flow problem matrix construction" begin
    bfp = test_bflow()
    (B,b) = BilevelFlowProblems.flow_constraint_standard(bfp)
    @test size(B)[1] == bfp.nv-2 + 1 + bfp.ne
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
    @test all(sum(B[i,:]) ≈ 1. for i in 4:8)
    @test all(b[4:8] .≈ (3., 2., 3., 2., 1.))
    @test B[4,5]  ≈ 1.
    @test B[5,9]  ≈ 1.
    @test B[6,13] ≈ 1.
    @test B[7,14] ≈ 1.
    @test B[8,15] ≈ 1.
end

@testset "Bilevel flow JuMP model" begin
    bfp = test_bflow()
    build_blp_model(bfp, CbcSolver())
end
