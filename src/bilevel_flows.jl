"""
Data structure and resolution methods of some forms of bilevel flow problems
"""
module BilevelFlowProblems

using SparseArrays: spzeros
using JuMP: @variable, @constraint, @objective
import JuMP

import ..MT
import ..VT

export BilevelFlowProblem

import ..build_blp_model
import ..SOS1Complementarity

"""
    BilevelFlowProblem{M1<:MT,M2<:AbstractMatrix{Bool},M3<:MT,M4<:AbstractArray{<:Real,3}}

Holding the information of a bilevel flow problem:
* Upper-level picks arcs to tax and a discrete tax level
* Lower-level solves a mincost flow problem with `minflow`
units out of the source and arc cost (init_cost + tax_options * chosen_option)
"""
struct BilevelFlowProblem{M1<:MT,M2<:AbstractMatrix{Bool},M3<:MT,M4<:AbstractArray{<:Real,3}}
    init_cost::M1
    taxable_edges::M2
    capacities::M3
    tax_options::M4
    nv::Int
    ne::Int
    minflow::Float64
    function BilevelFlowProblem(init_cost::M1,taxable_edges::M2,capacities::M3,tax_options::M4, minflow) where {M1<:MT,M2<:AbstractMatrix{Bool},M3<:MT,M4<:AbstractArray{<:Real,3}}
        size(init_cost) == size(taxable_edges) == size(capacities) == size(tax_options)[1:end-1] || DimensionMismatch("Matrix dimensions")
        nv = size(init_cost)[1]
        nv == size(init_cost)[2] || DimensionMismatch("Cost matrix should be square")
        ne = sum(1 for i in Base.OneTo(nv) for j in Base.OneTo(nv) if capacities[i,j] > 0. && i!=j)
        new{M1,M2,M3,M4}(init_cost, taxable_edges, capacities, tax_options, nv, ne, minflow)
    end
end

"""
	`build_blp_model(bfp::BilevelFlowProblem, solver; comp_method)`

Build the JuMP model from the data of `bfp` and assign it the passed solver.
Return `(m, r, y, f, λ, s)`
- `m`: `JuMP.Model`
- `y[i,j,k]`: binary variables indicating which option tax option is chosen by the upper-level
- `r[i,j]`: auxiliary variable storing the revenue made on arc `[i,j]`
- `f[i,j]`: flow variables (lower level)
- `λ`: flattened dual vector
- `s`: flattened slack variables for the flow problem
"""
function build_blp_model(bfp::BilevelFlowProblem, optimizer; comp_method = SOS1Complementarity())
    m = JuMP.Model(optimizer)
    (nv, _, nopt) = size(bfp.tax_options)
    @variable(m, y[i=1:nv,j=1:nv,k=1:nopt], Bin)
    @variable(m, r[i=1:nv,j=1:nv] >= 0)
    @variable(m, f[i=1:nv,j=1:nv] >= 0)
    @constraint(m, revenue_limit[i=1:nv,j=1:nv,k=1:nopt], r[i,j] <= bfp.tax_options[i,j,k] * f[i,j] + maximum(bfp.tax_options[i,j,:]) * bfp.capacities[i,j]*(1-y[i,j,k]))
    @constraint(m, unique_opt[i=1:nv,j=1:nv], sum(y[i,j,k] for k in 1:nopt) == 1.)
    @objective(m, Max,
        sum(r[i,j] for i in 1:nv for j in 1:nv if bfp.taxable_edges[i,j])
    )
    flat_flow = [f[i] for i in eachindex(f)]
    lin_cost  = [bfp.init_cost[i,j] + sum(y[i,j,k] * bfp.tax_options[i,j,k] for k in Base.OneTo(nopt)) for j in Base.OneTo(bfp.nv) for i in Base.OneTo(bfp.nv)]
    (B, b) = flow_constraint_standard(bfp)
    @variable(m,
        s[i=1:length(b)] >= 0
    )
    @constraint(m, B*flat_flow .+ s .== b)
    (_, λ, _) = build_blp_model(m, B, lin_cost, s, comp_method = comp_method)
    return (m, r, y, f, λ, s)
end

"""
    flow_constraint_standard(bfp::BilevelFlowProblem)

Construct the constraint matrix `B` and right-hand side vector `b`
for the lower-level problem
"""
function flow_constraint_standard(bfp::BilevelFlowProblem)
    B = spzeros( # constraint matrix
        1 +        # min flow from source
        bfp.nv-2 + # flow conservation at all nodes but source & sink
        bfp.nv*bfp.nv,    # capacity at all edges
        bfp.nv*bfp.nv
    )
    b = zeros(size(B)[1])
    b[1] = -bfp.minflow
    if bfp.capacities[1,bfp.nv] > 0.
        B[1,end-bfp.nv+1] = -1
    end
    for k in 2:bfp.nv-1
        if bfp.capacities[1,k] > 0.
            B[1,bfp.nv*(k-1)+1] = -1. # arc s -> k
        end
        for i in Base.OneTo(bfp.nv)
            if i != k
                if bfp.capacities[i,k] > 0.
                    B[k,bfp.nv*(k-1)+i] = 1.
                end
                if bfp.capacities[k,i] > 0.
                    B[k,bfp.nv*(i-1)+k] = -1.
                end
            end
        end
    end
    for i in Base.OneTo(bfp.nv), j in Base.OneTo(bfp.nv)
        B[bfp.nv-1+i+bfp.nv*(j-1),i+bfp.nv*(j-1)] = 1.
        b[bfp.nv-1+i+bfp.nv*(j-1)] = bfp.capacities[i,j]
    end
    return (B, b)
end

end
