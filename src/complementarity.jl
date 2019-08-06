"""
Method used for handling complementarity constraints `s[i] ⋅ λ[i] = 0`
Must implement:
```
add_complementarity_constraint(m, cm::ComplementarityMethod, s, λ)
```
"""
abstract type ComplementarityMethod end

"""
Handle complementarity constraints with
Special Ordered Sets of type 1
"""
struct SOS1Complementarity <: ComplementarityMethod end

struct BoundComplementarity{MD,MP} <: ComplementarityMethod
    Md::MD
    Mp::MP
    function BoundComplementarity(md::MD, mp::MP) where {MD <: Union{AbstractVector,Real}, MP <: Union{AbstractVector,Real}}
        return new{MD,MP}(md,mp)
    end
end

"""
Add complementarity constraints to the JuMP model `m` for each pair
`(s[i],λ[i])` with the corresponding complementarity method.
```
add_complementarity_constraint(m, cm::ComplementarityMethod, s, λ, y, σ)
```
"""
function add_complementarity_constraint end

"""
    add_complementarity_constraint(m, ::SOS1Complementarity, s, λ, y, σ)

Implements complementarity constraints using special ordered sets 1 (SOS1).
"""
function add_complementarity_constraint(m, ::SOS1Complementarity, s, λ, y, σ)
    @constraint(m, rhs_sos_constraints[i=eachindex(s)], [λ[i], s[i]] in MOI.SOS1([0.4, 0.6]))
    @constraint(m, nonnegative_sos_constraints[j=eachindex(y)], [σ[j], y[j]] in MOI.SOS1([0.4, 0.6]))
    return nothing
end

"""
Implements complementarity constraints using big-M type constraints
Introduces binary variables:
```
active_constraint[i] == 0 <=> λ[i] == 0
active_constraint[i] == 1 <=> s[i] == 0
variable_bound[j] == 0 <=> σ[j] == 0
variable_bound[j] == 1 <=> y[j] == 0
```
"""
function add_complementarity_constraint(m, bc::BoundComplementarity, s, λ, y, σ)
    ml = length(s)
    nl = length(y)
    @variable(m, active_constraint[1:ml], Bin)
    @variable(m, variable_bound[1:ml], Bin)
    add_bigm_primalbounds(m, bc, active_constraint, variable_bound, s, y)
    add_bigm_dualbounds(m, bc, active_constraint, variable_bound, λ, σ)
    return (active_constraint, variable_bound)
end

"""
Add primal bounds with one bound for each element
"""
function add_bigm_primalbounds(m, bc::BoundComplementarity{MD,MP}, active_constraint, variable_bound, s, y) where {MD,MP<:AbstractVector{<:Real}}
    ml = length(s)
    length(bc.Mp) == ml + length(y) || throw(DimensionMismatch("Primal bound vector"))
    @constraint(m, [i=1:ml],
        s[i] <= bc.Mp[i] * (1. - active_constraint[i])
    )
    @constraint(m, [j=eachindex(y)],
        y[j] <= bc.Mp[j+ml] * (1. - variable_bound[j])
    )
    return nothing
end

"""
Add primal bounds with the same bound for all elements
"""
function add_bigm_primalbounds(m, bc::BoundComplementarity{MD,MP}, active_constraint, variable_bound, s, y) where {MD,MP<:Real}
    @constraint(m, [i=eachindex(s)],
        s[i] <= bc.Mp * (1. - active_constraint[i])
    )
    @constraint(m, [j=eachindex(y)],
        y[j] <= bc.Mp * (1. - variable_bound[j])
    )
    return nothing
end

"""
Add dual bounds with one bound for each element
"""
function add_bigm_dualbounds(m, bc::BoundComplementarity{MD,<:Any}, active_constraint, variable_bound, λ, σ) where {MD<:AbstractVector{<:Real}}
    ml = length(λ)
    length(bc.Md) == ml + length(σ) || throw(DimensionMismatch("Dual bound vector"))
    @constraint(m, [i=1:ml],
        λ[i] <= bc.Md[i] * active_constraint[i]
    )
    @constraint(m, [j=1:length(σ)],
        σ[j] <= bc.Md[ml+j] * variable_bound[j]
    )
    return nothing
end

"""
Add dual bounds with the same bound for all elements
"""
function add_bigm_dualbounds(m, bc::BoundComplementarity{MD,<:Any}, active_constraint, variable_bound, λ, σ) where {MD<:Real}
    @constraint(m, [i=1:length(λ)],
        λ[i] <= bc.Md * active_constraint[i]
    )
    @constraint(m, [j=1:length(σ)],
        σ[j] <= bc.Md * variable_bound[j]
    )
    return nothing
end
