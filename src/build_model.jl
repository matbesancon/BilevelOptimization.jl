
"""
Build a JuMP model (constraints and objective)
based on the data from bp
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

function build_blp_model(m::JuMP.Model, B::M, d, x, y, s, F) where {M<:MT}
    ml = size(B)[1]
    @variable(m, λ[1:ml] >= 0)
    if F !== nothing
        @constraint(m, d .+ F' * x .+ B' * λ .== 0.0)
    else
        @constraint(m, d .+ B' * λ .== 0.0)
    end

    for i in Base.OneTo(ml)
        JuMP.addSOS1(m, [λ[i], s[i]])
    end
    return (m, λ)
end

function build_blp_model(m::JuMP.Model, B::M, d, y, s) where {M<:MT}
    ml = size(B)[1]
    @variable(m, λ[1:ml] >= 0)
    @constraint(m, d .+ B' * λ .== 0.0)
    for i in Base.OneTo(ml)
        JuMP.addSOS1(m, [λ[i], s[i]])
    end
    return (m, λ)
end
