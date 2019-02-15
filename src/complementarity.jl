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
        new{MD,MP}(Md,Mp)
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
Implements complementarity constraints using special ordered sets 1
Called with either
```
add_complementarity_constraint(m, SOS1Complementarity, s, λ, y, σ)
add_complementarity_constraint(m, SOS1Complementarity(), s, λ, y, σ)
```
"""
function add_complementarity_constraint(m, ::S, s, λ, y, σ) where {S <: Union{SOS1Complementarity,Type{SOS1Complementarity}}}
    for i in eachindex(s)
        JuMP.addSOS1(m, [λ[i], s[i]])
    end
    for j in eachindex(y)
        JuMP.addSOS1(m, [σ[j], y[j]])
    end
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
    length(bc.mp) == ml + length(y) || throw(DimensionMismatch("Primal bound vector"))
    @constraint(m, [i=1:ml],
        s[i] <= bc.mp[i] * (1. - active_constraint[i])
    )
    @constraint(m, [j=1:length(y)],
        y[j] <= bc.mp[j+ml] * (1. - variable_bound[j])
    )
    return nothing
end

"""
Add primal bounds with the same bound for all elements
"""
function add_bigm_primalbounds(m, bc::BoundComplementarity{MD,MP}, active_constraint, variable_bound, s, y) where {MD,MP<:Real}
    @constraint(m, [i=1:length(s)],
        s[i] <= bc.mp * (1. - active_constraint[i])
    )
    @constraint(m, [j=1:length(y)],
        y[j] <= bc.mp * (1. - variable_bound[j])
    )
    return nothing
end

"""
Add dual bounds with one bound for each element
"""
function add_bigm_dualbounds(m, bc::BoundComplementarity{MD,MP}, active_constraint, variable_bound, λ, σ) where {MD<:AbstractVector{<:Real},MP}
    ml = length(λ)
    length(bc.md) == ml + length(σ) || throw(DimensionMismatch("Dual bound vector"))
    @constraint(m, [i=1:ml],
        λ[i] <= bc.md[i] * active_constraint[i]
    )
    @constraint(m, [j=1:length(σ)],
        σ[j] <= bc.md[ml+j] * variable_bound[j]
    )
    return nothing
end

"""
Add dual bounds with the same bound for all elements
"""
function add_bigm_dualbounds(m, bc::BoundComplementarity{MD,MP}, active_constraint, variable_bound, λ) where {MD<:Real,MP}
    @constraint(m, [i=1:length(λ)],
        λ[i] <= bc.md * active_constraint[i]
    )
    @constraint(m, [j=1:length(σ)],
        σ[j] <= bc.md * variable_bound[j]
    )
    return nothing
end
