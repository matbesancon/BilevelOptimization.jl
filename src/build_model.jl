
"""
Build a JuMP model (constraints and objective)
based on the data from bp and with a solver
"""
function build_blp_model(bp::BilevelLP, solver)
    m = JuMP.Model(solver = solver)
    @variable(m, x[1:bp.nu] >= 0)
    @variable(m, y[1:bp.nl] >= 0)
    @constraint(m, uppercons[i=1:bp.mu],
        sum(bp.G[i,j]*x[j] for j in 1:bp.nu) +
        sum(bp.H[i,j]*y[j] for j in 1:bp.nl) <= bp.q[i])
    @objective(m, Min, sum(bp.cx .* x) + sum(bp.cy .* y))
    for j in bp.Jx
        JuMP.setcategory(x[j], :Int)
    end
    # adding SOS1 constraints
    return build_blp_model(m, bp, x, y)
end

"""
Add the lower-level constraints and optimality conditions to
an existing JuMP model. This assumes the upper-level feasibility
constraints and objective have already been set
"""
function build_blp_model(m::JuMP.Model, bp::BilevelLP, x, y)
    @variable(m, s[1:bp.ml] >= 0) # lower-level slack variables
    @constraint(m, lowercons[i=1:bp.ml],
        sum(bp.A[i,j]*x[j] for j in 1:bp.nu) +
        sum(bp.B[i,j]*y[j] for j in 1:bp.nl) + s[i] == bp.b[i]
    )
    @variable(m, λ[1:bp.ml] >= 0)
    @constraint(m, bp.d .+ bp.F' * x .+ bp.B' * λ .== 0.0)
    for i in Base.OneTo(bp.ml)
        JuMP.addSOS1(m, [λ[i], s[i]])
    end
    return (m, x, y, λ)
end

"""
Build the bilevel JuMP model from the data without
grouping everything into a `BilevelLP`
"""
function build_blp_model(m::JuMP.Model, B::M, d, x, y, s, F) where {M<:MT}
    ml = size(B)[1]
    @variable(m, λ[1:ml] >= 0)
    @constraint(m, d .+ F' * x .+ B' * λ .== 0.0)
    for i in Base.OneTo(ml)
        JuMP.addSOS1(m, [λ[i], s[i]])
    end
    return (m, λ)
end

function build_blp_model(m::JuMP.Model, B::M, d, s) where {M<:MT}
    (ml,nl) = size(B)
    @variable(m, λ[1:ml] >= 0)
    @constraint(m, kkt, d .+ B' * λ .== 0.0)
    for i in Base.OneTo(ml)
        JuMP.addSOS1(m, [λ[i], s[i]])
    end
    return (m, λ)
end
