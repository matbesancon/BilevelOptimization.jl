
"""
Convenience alias for a matrix of a sub-type of `Real`
"""
const MT = AbstractMatrix{<:Real}
"""
Convenience alias for a vector of a sub-type of `Real`
"""
const VT = AbstractVector{<:Real}

"""
A bilevel linear optimization problem of the form:
```
min cx^T * x + cy^T * y
s.t. G x + H y <= q
     x_j ∈ [xl_j,xu_j]
     x_j ∈ ℤ ∀ j ∈ Jx
     y ∈ arg min {
        d^T * y + x^T * F * y
        s.t. A x + B y <= b
             y_j ∈ [yl_j,yu_j]
        }
```
Note that integer variables are allowed at the upper level.
"""
mutable struct BilevelLP{V<:VT,M<:MT,MQ<:MT,VB<:AbstractVector{Bool}}
    cx::V
    cy::V
    G::M
    H::M
    q::V
    d::V
    A::M
    B::M
    b::V
    nu::Int
    nl::Int
    mu::Int
    ml::Int
    xl::V
    xu::V
    yl::VB # if lowerbound yi >= 0, otherwise yi >= -∞. For other bounds, add a constraint
    yu::V
    Jx::Vector{Int} # ∀ j ∈ Jx, x[j] is integer
    F::MQ

    function BilevelLP(cx::V,
                       cy::V,
                       G::M,
                       H::M,
                       q::V,
                       d::V,
                       A::M,
                       B::M,
                       b::V,
                       Jx::Vector{Int} = Int[],
                       F::MQ = zeros(length(cx),length(cy)); ylowerbound = true) where {V<:VT,M<:MT,MQ<:MT}
        nu = length(cx)
        nl = length(cy)
        nl == length(d) || DimensionMismatch("Objectives")
        mu = length(q)
        ml = length(b)
        size(A) == (ml,nu) && size(B) == (ml,nl) || DimensionMismatch("Lower constraints")
        size(G) == (mu,nu) && size(H) == (mu,nl) || DimensionMismatch("Higher constraints")
        xl = zeros(nu)
        xu = Inf64 .* ones(nu)
        yl = zeros(Bool, nl) .+ ylowerbound # by default, y >= 0
        yu = Inf64 .* ones(nl)
        size(F) == (nu,nl) || DimensionMismatch("Quadratic constraint")
        new{V,M,MQ,Vector{Bool}}(cx,cy,G,H,q,d,A,B,b,nu,nl,mu,ml,xl,xu,yl,yu,Jx,F)
    end
end

"""
VariableType enum distinguishing upper- and lower-level
variables for setting upper and lower bounds
"""
@enum VariableType begin
    lower
    upper
end

"""
Set a lower bound on a lower or higher variable of `bp` depending on `vartype`
"""
function JuMP.setlowerbound(bp::BilevelLP, vartype::VariableType, j::Integer, v::T) where {T<:Real}
    if vartype == lower::VariableType
        if v ≈ 0.
            bp.yl[j] = true
        elseif isinf(v) && v < 0.
            bp.yl[j] = false
        else
            throw(DomainError(v,"The lower bound on a lower variable can only be 0 or -∞"))
        end
    else
        bp.xl[j] = v
    end
    return nothing
end

"""
Set an upper bound on a lower or higher variable of `bp` depending on `vartype`
Since we need the dual for lower-level upper bounds, this case is registered as
a constraint.
"""
function JuMP.setupperbound(bp::BilevelLP, vartype::VariableType, j::Integer, v::T) where {T<:Real}
    if vartype == upper::VariableType
        bp.xu[j] = v
    else
        row = zeros(1,bp.nl); row[1,j] = 1.
        bp.B = vcat(bp.B, row)
        bp.A = vcat(bp.A, zeros(1,bp.nu))
        push!(bp.b, v)
        bp.ml = bp.ml + 1
    end
end
