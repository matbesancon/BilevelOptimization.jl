using BilevelOptimization
using Test

import Cbc
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
    (m, x, y, λ) = build_blp_model(bp, Cbc.CbcSolver())
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
    (m, x, y, λ) = build_blp_model(bp, Cbc.CbcSolver())
    @test JuMP.getcategory(x[2]) == :Int
    @test JuMP.getcategory(x[1]) == :Cont
    @test JuMP.getcategory(y[1]) == :Cont
end
