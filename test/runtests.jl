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
    (m, x, y, f, λ) = build_blp_model(bfp, CbcSolver())
    st = JuMP.solve(m)
    @test st === :Optimal
    @test getobjectivevalue(m) ≈ 6.
    for j in 1:size(x)[2]
        for i in 1:size(x)[1]
            @test getvalue(x[i,j]) ≈ sum(getvalue(y[i,j,:]).*bfp.tax_options[i,j,:]) * getvalue(f[i,j])
        end
    end
end
